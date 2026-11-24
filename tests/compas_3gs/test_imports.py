import os
import subprocess
import sys

SRC = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", "src"))


def test_import_utilities_without_matplotlib():
    code = "import sys; import compas_3gs.utilities; assert 'matplotlib' not in sys.modules, 'matplotlib was imported'"
    env = dict(os.environ, PYTHONPATH=SRC)
    result = subprocess.run([sys.executable, "-c", code], env=env, capture_output=True, text=True, check=False)
    assert result.returncode == 0, result.stderr
