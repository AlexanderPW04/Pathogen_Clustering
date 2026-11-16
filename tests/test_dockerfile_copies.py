"""The Docker image must contain every top-level module the web app imports.

The Dockerfile copies top-level files one by one, so a new module (such as
results_layout.py) is easy to forget; the hosted worker then fails to boot
with ModuleNotFoundError (found on Render, 2026-11-16).
"""

import re
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

import app  # noqa: E402,F401  (imports everything the web app needs)


def copied_by_dockerfile():
    names = set()
    for line in (ROOT / "Dockerfile").read_text(encoding="utf-8").splitlines():
        if line.startswith("COPY "):
            parts = line.split()[1:-1]  # the last part is the destination
            names.update(part.rstrip("/") for part in parts)
    return names


class DockerfileTest(unittest.TestCase):
    def test_every_module_the_app_imports_is_copied_into_the_image(self):
        copied = copied_by_dockerfile()
        needed = set()
        for module in list(sys.modules.values()):
            path = getattr(module, "__file__", None)
            if not path:
                continue
            path = Path(path).resolve()
            if ROOT not in path.parents or ".venv" in path.parts or "tests" in path.relative_to(ROOT).parts:
                continue
            relative = path.relative_to(ROOT)
            needed.add(relative.parts[0] if len(relative.parts) > 1 else relative.name)
        missing = sorted(name for name in needed if name not in copied)
        self.assertEqual(missing, [], f"add these to a COPY line in the Dockerfile: {missing}")
        self.assertIn("results_layout.py", copied)
        self.assertTrue(re.search(r"COPY genomics/ ./genomics/", (ROOT / "Dockerfile").read_text()))


if __name__ == "__main__":
    unittest.main()
