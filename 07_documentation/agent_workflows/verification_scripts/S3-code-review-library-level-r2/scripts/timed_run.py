# Runs one command from the current directory; reports wall time, CPU time and peak resident memory of the child,
# and samples the load (number of other R processes) once per second.
import resource, subprocess, sys, time, os
out_dir = sys.argv[1]; log = sys.argv[2]; cmd = sys.argv[3:]
env = dict(os.environ); env["S3_OUT_DIR"] = out_dir
t0 = time.time()
others = []
with open(log, "w") as fh:
    p = subprocess.Popen(cmd, stdout=fh, stderr=subprocess.STDOUT, env=env)
    while p.poll() is None:
        try:
            n = int(subprocess.run("pgrep -c -x R", shell=True, capture_output=True, text=True).stdout.strip() or 0)
        except Exception:
            n = -1
        others.append(n)
        time.sleep(1)
    rc = p.returncode
t1 = time.time()
ru = resource.getrusage(resource.RUSAGE_CHILDREN)
with open(log + ".timing", "w") as fh:
    fh.write("exit_code %d\nwall_seconds %.1f\nuser_seconds %.1f\nsys_seconds %.1f\nmax_rss_MB %.0f\nR_processes_seen_min %d\nR_processes_seen_max %d\n" % (
        rc, t1 - t0, ru.ru_utime, ru.ru_stime, ru.ru_maxrss / 1024.0, min(others) if others else -1, max(others) if others else -1))
