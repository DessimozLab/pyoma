import hashlib
import os
import logging

logger = logging.getLogger(__name__)


class JobArray(object):
    def __init__(self, nr_procs, this_proc_nr):
        self.nr_procs = int(nr_procs)
        self.this_proc_nr = int(this_proc_nr)
        if (
            self.nr_procs < 1
            or self.this_proc_nr < 1
            or self.this_proc_nr > self.nr_procs
        ):
            raise ValueError(
                "cannot determine HPC parallel job ids: nr_procs={}, this_proc_nr={}".format(
                    nr_procs, this_proc_nr
                )
            )

    def is_my_job(self, chunk):
        chunk = str(chunk).encode("utf-8") if not isinstance(chunk, bytes) else chunk
        h = hashlib.md5(chunk)
        as_int = int(h.hexdigest(), 16)
        res = (as_int % self.nr_procs) == (self.this_proc_nr - 1)
        logger.debug("chunk {} to be processed: {}".format(chunk, res))
        return res

    def __repr__(self):
        return "{}({},{})".format(
            self.__class__.__name__, self.nr_procs, self.this_proc_nr
        )

    def __str__(self):
        return "Jobarray process {} of {}".format(self.this_proc_nr, self.nr_procs)

    def modify_filename(self, basename):
        if self.nr_procs > 1:
            base, ext = os.path.splitext(basename)
            return "{}_{}-{}{}".format(base, self.this_proc_nr, self.nr_procs, ext)
        return basename


def detect_hpc_jobarray(nr_procs, this_proc_nr=None):
    if this_proc_nr is not None:
        return JobArray(nr_procs, int(this_proc_nr))
    if os.getenv("LSB_JOBID"):
        return JobArray(nr_procs, os.getenv("LSB_JOBINDEX", "1"))
    elif os.getenv("SGE_TASK_ID"):
        return JobArray(nr_procs, os.getenv("SGE_TASK_ID"))
    elif os.getenv("SLURM_ARRAY_JOB_ID"):
        return JobArray(nr_procs, os.getenv("SLURM_ARRAY_TASK_ID"))
    elif os.getenv("THIS_PROC_NR"):
        return JobArray(nr_procs, os.getenv("THIS_PROC_NR"))
    elif nr_procs is None or nr_procs == 1:
        return JobArray(1, 1)
    else:
        raise ValueError("Cannot identify HPC jobarray setup")
