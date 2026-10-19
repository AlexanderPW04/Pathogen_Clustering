import os
import subprocess
import unittest
from unittest.mock import Mock, patch

import requests

from geocode_locations import check_geonames
from system_diagnostics import check_mafft


class SystemDiagnosticsTests(unittest.TestCase):
    def test_mafft_success_aligns_two_sequences_and_reports_version(self):
        def fake_run(command, **kwargs):
            if "--version" in command:
                return subprocess.CompletedProcess(command, 0, stdout="", stderr="v7.505")
            return subprocess.CompletedProcess(
                command,
                0,
                stdout=">DIAG_A\nACGTACGTACGTACGTACGT\n"
                ">DIAG_B\nACGTACGTACG-ACGTACGT\n",
                stderr="",
            )

        with patch("system_diagnostics.find_mafft", return_value="/usr/bin/mafft"):
            with patch("system_diagnostics.subprocess.run", side_effect=fake_run) as run:
                result = check_mafft()

        self.assertEqual(result["status"], "Operational")
        self.assertEqual(result["version"], "v7.505")
        self.assertEqual(run.call_count, 2)
        self.assertEqual(run.call_args_list[0].args[0][:3], ["/usr/bin/mafft", "--auto", "--nuc"])
        self.assertEqual(run.call_args_list[0].kwargs["timeout"], 5)

    def test_mafft_missing_executable(self):
        with patch("system_diagnostics.find_mafft", side_effect=RuntimeError("missing")):
            result = check_mafft()
        self.assertEqual(result["status"], "Unavailable")

    def test_mafft_failed_alignment_and_timeout(self):
        with patch("system_diagnostics.find_mafft", return_value="/usr/bin/mafft"):
            with patch(
                "system_diagnostics.subprocess.run",
                return_value=subprocess.CompletedProcess("mafft", 1, stdout=""),
            ):
                failed = check_mafft()
            with patch(
                "system_diagnostics.subprocess.run",
                side_effect=subprocess.TimeoutExpired("mafft", 5),
            ):
                timed_out = check_mafft()
        self.assertEqual(failed["status"], "Check failed")
        self.assertEqual(timed_out["status"], "Check failed")

    def test_mafft_rejects_missing_or_unequal_sequences(self):
        with patch("system_diagnostics.find_mafft", return_value="/usr/bin/mafft"):
            for output in (">DIAG_A\nACGT\n", ">DIAG_A\nACGT\n>DIAG_B\nACG\n"):
                with self.subTest(output=output):
                    with patch(
                        "system_diagnostics.subprocess.run",
                        return_value=subprocess.CompletedProcess("mafft", 0, stdout=output),
                    ):
                        result = check_mafft()
                    self.assertEqual(result["status"], "Check failed")

    def test_mafft_windows_batch_file_uses_cmd(self):
        def fake_run(command, **kwargs):
            if "--version" in command:
                return subprocess.CompletedProcess(command, 1, stdout="")
            return subprocess.CompletedProcess(
                command, 0,
                stdout=">DIAG_A\nACGTACGTACGTACGTACGT\n"
                ">DIAG_B\nACGTACGTACG-ACGTACGT\n",
            )

        with patch("system_diagnostics.find_mafft", return_value=r"C:\MAFFT\mafft.bat"):
            with patch("system_diagnostics.subprocess.run", side_effect=fake_run) as run:
                result = check_mafft()
        self.assertEqual(result["status"], "Operational")
        self.assertEqual(run.call_args_list[0].args[0][:2], ["cmd.exe", "/c"])

    def test_geonames_missing_configuration_makes_no_request(self):
        with patch.dict(os.environ, {"GEONAMES_USERNAME": ""}):
            with patch("geocode_locations.requests.get") as get:
                result = check_geonames()
        self.assertEqual(result["status"], "Not configured")
        get.assert_not_called()

    def test_geonames_success_uses_secure_endpoint_and_valid_coordinates(self):
        response = Mock()
        response.json.return_value = {
            "geonames": [{"lat": "33.4255", "lng": "-111.9400"}]
        }
        with patch.dict(os.environ, {"GEONAMES_USERNAME": "private-test-account"}):
            with patch("geocode_locations.requests.get", return_value=response) as get:
                result = check_geonames()

        self.assertEqual(result["status"], "Operational")
        self.assertNotIn("private-test-account", str(result))
        get.assert_called_once()
        self.assertEqual(get.call_args.args[0], "https://secure.geonames.org/searchJSON")
        self.assertIs(get.call_args.kwargs["verify"], True)
        self.assertEqual(get.call_args.kwargs["timeout"], 4)

    def test_geonames_connection_and_api_errors_are_distinct(self):
        response = Mock()
        response.json.return_value = {"status": {"message": "private API detail"}}
        with patch.dict(os.environ, {"GEONAMES_USERNAME": "private-test-account"}):
            with patch(
                "geocode_locations.requests.get", side_effect=requests.ConnectionError()
            ):
                connection = check_geonames()
            with patch("geocode_locations.requests.get", return_value=response):
                api_error = check_geonames()

        self.assertEqual(connection["status"], "Unavailable")
        self.assertEqual(api_error["status"], "Check failed")
        self.assertNotIn("private API detail", str(api_error))

    def test_geonames_http_api_error(self):
        response = Mock()
        response.raise_for_status.side_effect = requests.HTTPError("private detail")
        with patch.dict(os.environ, {"GEONAMES_USERNAME": "private-test-account"}):
            with patch("geocode_locations.requests.get", return_value=response):
                result = check_geonames()
        self.assertEqual(result["status"], "Check failed")
        self.assertNotIn("private detail", str(result))

    def test_geonames_rejects_out_of_range_coordinates(self):
        response = Mock()
        response.json.return_value = {"geonames": [{"lat": "200", "lng": "-111.94"}]}
        with patch.dict(os.environ, {"GEONAMES_USERNAME": "private-test-account"}):
            with patch("geocode_locations.requests.get", return_value=response):
                result = check_geonames()
        self.assertEqual(result["status"], "Check failed")


if __name__ == "__main__":
    unittest.main()
