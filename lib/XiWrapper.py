import os
import subprocess
import sys  # get the number of command line arguments
import argparse     # parse command line arguments
import logging
import time
import datetime


logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())


class XiSearchException(subprocess.CalledProcessError):
    def __init__(self, returncode, cmd, out_file, output=None):
        subprocess.CalledProcessError.__init__(self, returncode, cmd, output)
        self.out_file = out_file
        pass

    def __str__(self):
        return "XiSearch command '{}' produced an error: '{}'"\
            .format(" ".join(self.cmd), self.output)


class XiSearchOutOfMemoryException(XiSearchException):
    def __init__(self, returncode, cmd, out_file, output):
        XiSearchException.__init__(self, returncode, cmd, out_file, output)
        pass

    def __str__(self):
        return "Command '{:s}' produced java Memory Exception '{}'. You might want to remove the result stub."\
            .format(" ".join(self.cmd), self.output)


class XiSearchDaemoniseFailureException(XiSearchException):
    def __init__(self, returncode, cmd, out_file, output):
        XiSearchException.__init__(self, returncode, cmd, out_file, output)
        pass

    def __str__(self):
        return "Command '{:s}' failed while writing result with following error: '{}'." + \
               "The result file should be mostly complete though."\
            .format(" ".join(self.cmd), self.output)


class XiWrapper:
    """
    XiWrapper class to run Xi searches from command line.
    """

    def __init__(self):
        pass

    @staticmethod
    def calculate_elapsed_time(starttime):
        """ returns elapsed time since starttime"""
        seconds_passed = (time.time() - starttime)
        return str(datetime.timedelta(seconds=seconds_passed))

    @staticmethod
    def build_xi_arguments(xi_path, xi_config, peak_files, fasta_files, memory, output, additional_parameters=()):
        """
        Example command (fastutil-version for mem-optimization has to match jave version):
        java -cp XiSearch.jar:fastutil-8.1.0.jar rappsilber.applications.Xi --config=[path to config] --xiconf=TOPMATCHESONLY:true --peaks=[path to peaklist1] --peaks=[path to peaklist2] --peaks=[path to peaklist3] --fasta=[path to fasta file1] --fasta=[path to fasta file2] --fasta=[path to fasta file3] --output=[path to result file]
        """
        cmd = ["java"]

        # memory
        if memory:
            cmd.extend(["-Xmx" + memory
                        # , "-XX:NativeMemoryTracking=detail"    # for java memory tracking
                        ])

        # XiSearch jar
        cmd.extend([
            "-cp"
            , xi_path
            # Win:
            # + ";fastutil-7.2.1.jar"    # uncomment for java 7 mem optimization
            # + ";fastutil-8.1.0.jar"    # uncomment for java 8 mem optimization
            # linux:
            + ":fastutil-7.2.1.jar"    # uncomment for java 7 mem optimization
            # + ":fastutil-8.1.0.jar"    # uncomment for java 8 mem optimization
            , "rappsilber.applications.Xi"
        ])

        # config
        cmd.append("--config=" + xi_config)

        # additional parameters
        for par in additional_parameters:
            cmd.append(par)

        # peak_files
        for peak_file in peak_files:
            cmd.append("--peaks=" + peak_file)

        # fasta_files
        for fasta_file in fasta_files:
            cmd.append("--fasta=" + fasta_file)

        # output file
        cmd.append("--output=" + output)
        return cmd

    @staticmethod
    def xi_execution(xi_config, peak_files, fasta_files, memory=None, output_file="xi_results",
                     additional_parameters=list(),
                     xi_path="XiSearch.jar"):
        """
        Calls Xi and gives back Xi result filepath
        Xi gives back a single csv file

        Example cmd line XiSearch execution:
        # java -cp XiSearch.jar rappsilber.applications.Xi --conf=[path to config]
        # --xiconfig=TOPMATCHESONLY:true --peaks=[path to peaklist1] --peaks=[path to peaklist2] --peaks=[path to peaklist3]
        # --fasta=[path to fasta file1] --fasta=[path to fasta file2] --fasta=[path to fasta file3]
        # --output=[path to result file]

        :param xi_config: string
        :param peak_files: list of strings
        :param fasta_files: list of strings
        :param memory: string of format int[gmk], i.e. "30g" for 30GB of RAM
        :param output_file: string
        :param additional_parameters: list of strings, optional
        :param xi_path: path to xisearch jar file, optional
        :return: string: output_file.csv
        """
        assert type(peak_files) == type(fasta_files) == type(additional_parameters) == list, \
            """type of the following files needs to be list. It is actually:
            peak_files: {}
            fasta_files: {}
            additional_parameters: {}""" \
                .format(type(peak_files), type(fasta_files), type(additional_parameters))

        list_of_all_files = peak_files + fasta_files
        list_of_all_files.append(xi_path)
        for f in list_of_all_files:
            if not os.path.exists(f):
                raise IOError("Could not find specified file: '{}' Is the path correct?"
                              .format(f))

        if not os.path.exists(os.path.split(output_file)[0]):
            os.makedirs(os.path.split(output_file)[0])
        output_file = os.path.splitext(output_file)[0] + ".csv"

        # generate xi commands
        xi_cmd = XiWrapper.build_xi_arguments(xi_path=xi_path,
                                              xi_config=xi_config,
                                              peak_files=peak_files,
                                              fasta_files=fasta_files,
                                              memory=memory,
                                              output=output_file,
                                              additional_parameters=additional_parameters)

        # call xi
        starttime = time.time()
        logger.info("XiSearch cmd: {}".format(" ".join(map(str, xi_cmd))))
        process = subprocess.Popen(xi_cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        # real time output of Xi messages
        while True:
            output = process.stdout.readline()
            exit_code = process.poll()
            if output == '' and exit_code is not None:
                break
            elif output:
                # print output.strip()
                logger.debug("XiSearch: " + output.strip())
                if "java.lang.OutOfMemoryError" in output:
                    process.kill()
                    raise XiSearchOutOfMemoryException(returncode=1, cmd=xi_cmd, out_file=output_file, output=output)
                elif "could not daemonise BufferedResultWriter_batchforward" in output:
                    process.kill()
                    raise XiSearchDaemoniseFailureException(
                        returncode=1, cmd=xi_cmd, out_file=output_file, output=output
                    )

        if exit_code != 0:  # if process exit code is non zero
            raise XiSearchException(exit_code, xi_cmd, output_file, 'XiSearch exited with error message!')
        logger.info("XiSearch execution took {} for cmd: {}"
                    .format(XiWrapper.calculate_elapsed_time(starttime), xi_cmd))
        return output_file
