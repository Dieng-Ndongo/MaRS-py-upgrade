import os
import sys
import subprocess
from pathlib import Path

BASE_DIR = Path(__file__).resolve().parent
PIPELINE_SCRIPT = BASE_DIR / "pipeline_python.py"

ENV_NAME = "pipeline_env"

def run_cmd_live(cmd, env=None):
    print(f"[CMD] {' '.join(cmd)}")

    process = subprocess.Popen(
        cmd,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        env=env,
        cwd=str(BASE_DIR)
    )

    for line in process.stdout:
        print(line, end="")

    process.wait()

    if process.returncode != 0:
        raise subprocess.CalledProcessError(process.returncode, cmd)

def main():
    print("Step 1/3: Préparation environnement...")

    # 👉 Si tu veux Docker, remplace ici
    # cmd = ["docker", "run", ...]
    
    print("Step 2/3: Lancement pipeline...")

    cmd = ["python", str(PIPELINE_SCRIPT)]
    run_cmd_live(cmd)

    print("Step 3/3: Terminé ✔")

if __name__ == "__main__":
    main()