import os
import subprocess
import time
import logging
import datetime


logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())


def calculate_elapsed_time(starttime):
    """ returns elapsed time since starttime in minutes """
    seconds = (time.time() - starttime)
    return str(datetime.timedelta(seconds=seconds))


class XiFdrWrapper:
    def __init__(self):
        pass

    @staticmethod
    def build_xifdr_arguments(fdr_input_csv, fdr_output_dir, pepfdr, memory="1G", reportfactor="10000",
                              additional_xifdr_arguments=list(), xifdr_filename="xiFDRDB-1.1.25.55-jar-with-dependencies.jar"):
        assert isinstance(fdr_input_csv, (list, tuple, str, unicode)), \
            """type of fdr_input_csv needs to be in (list, tuple, str, unicode) but is: {}""".format(type(fdr_input_csv))
        if isinstance(fdr_input_csv, (str, unicode)):
            fdr_input_csv = [fdr_input_csv]
        # Example cmd line:
        # java -Xmx1g -cp xiFDRDB-1.0.13.32-jar-with-dependencies.jar org.rappsilber.fdr.CSVinFDR --psmfdr=X --pepfdr=X
        # --proteinfdr=X --reportfactor=X --linkfdr=X --ppifdr=X --csvOutDir=X --csvBaseName=X csv-file1 csv-file2
        cmd = []
        # memory
        cmd.extend(["java", "-Xmx" + memory, "-cp", xifdr_filename,
                    "org.rappsilber.fdr.CSVinFDR", '--reportfactor='+reportfactor])
        # pepfdr
        cmd.append("--pepfdr=" + pepfdr)
        # additional arguments
        for par in additional_xifdr_arguments:
            cmd.append(par)
        # reportfactor
        # taken into args as it will probably never be changed
        # cmd.append("--reportfactor=" + reportfactor)
        # csvOutDir
        cmd.append("--csvOutDir=" + fdr_output_dir)
        # input csv
        for i in fdr_input_csv:
            cmd.append(i)
        return cmd

    @staticmethod
    def xifdr_execution(
            xifdr_input_csv, xifdr_output_dir, pepfdr="5", memory="1G", reportfactor="10000",
            additional_xifdr_arguments=list(), xifdr_filename="xiFDRDB-1.1.25.55-jar-with-dependencies.jar"
    ):
        """
        takes XiSearch output and gives back fdr csv
        :param xifdr_input_csv:
        :param xifdr_output_dir:
        :param additional_xifdr_arguments:
        :return:
        list of result files from xifdr
        """
        list_of_results = []
        if isinstance(xifdr_input_csv, str):
            xifdr_input_csv = [xifdr_input_csv]
        if not os.path.exists(xifdr_output_dir):
            os.makedirs(xifdr_output_dir)
        starttime = time.time()
        xifdr_cmd = XiFdrWrapper.build_xifdr_arguments(xifdr_input_csv, xifdr_output_dir, pepfdr, memory, reportfactor,
                                                       additional_xifdr_arguments, xifdr_filename=xifdr_filename)
        logger.info("xiFDR arguments: {}".format(" ".join(map(str, xifdr_cmd))))
        # # # # #
        process = subprocess.Popen(xifdr_cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        # real time output of Xi messages
        while True:
            output = process.stdout.readline()
            exit_code = process.poll()
            if output == '' and exit_code is not None:
                break
            if output:
                # print output.strip()
                logger.debug("xiFDR: " + output.strip())
        if exit_code != 0:  # if process exit code is non zero
            raise subprocess.CalledProcessError(exit_code, " ".join(xifdr_cmd))
        logger.info("xiFDR execution took {} for cmd: {}"
                      .format(calculate_elapsed_time(starttime), xifdr_cmd))
        # read filenames from result dir
        for rel_dir, sub_dirs, files in os.walk(xifdr_output_dir):
            list_of_results = [os.path.join(rel_dir, f) for f in files]
        return list_of_results
