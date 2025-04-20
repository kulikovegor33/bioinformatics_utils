import pytest
import os
import logging
from bioinformatics_utils import filter_fastq
from pathlib import Path

class TestFileIO:

    def test_filter_fastq_creates_output(self, tmp_path):
        input_file = tmp_path / "input.fastq"
        output_file = tmp_path / "output.fastq"

        input_file.write_text("""@seq1
ATGCATGC
+
IIIIIIII
""")

        filter_fastq(str(input_file), str(output_file))
        assert output_file.exists()
        assert output_file.read_text().strip() != ""

    def test_empty_output_if_all_filtered(self, tmp_path):
        input_file = tmp_path / "low_quality.fastq"
        output_file = tmp_path / "filtered.fastq"

        input_file.write_text("""@seq1
ATGCATGC
+
########
""")  # low quality ASCII 35

        filter_fastq(str(input_file), str(output_file), quality_threshold=50)
        assert output_file.read_text().strip() == ""


class TestFilteringLogic:

    def test_gc_content_filtering(self, tmp_path):
        input_file = tmp_path / "gc_test.fastq"
        output_file = tmp_path / "gc_out.fastq"

        input_file.write_text("""@seq1
ATATATAT
+
IIIIIIII
""")  # GC = 0%

        filter_fastq(str(input_file), str(output_file), gc_bounds=(50, 100))
        assert output_file.read_text().strip() == ""

    def test_length_filtering(self, tmp_path):
        input_file = tmp_path / "length.fastq"
        output_file = tmp_path / "length_out.fastq"

        input_file.write_text("""@seq1
ATGCATGCATGC
+
IIIIIIIIIIII
""")

        filter_fastq(str(input_file), str(output_file), length_bounds=(0, 5))
        assert output_file.read_text().strip() == ""


class TestErrorHandling:

    def test_input_file_not_found(self):
        with pytest.raises(FileNotFoundError):
            filter_fastq("non_existent_file.fastq", "out.fastq")

    def test_invalid_format(self, tmp_path):
        input_file = tmp_path / "bad.fastq"
        output_file = tmp_path / "bad_out.fastq"

        input_file.write_text("This is not FASTQ at all")

        with pytest.raises(Exception):
            filter_fastq(str(input_file), str(output_file))


class TestLogFile:

    @pytest.fixture(autouse=True)
    def cleanup_log(self, tmp_path):
        self.log_path = tmp_path / "filter.log"
        if self.log_path.exists():
            self.log_path.unlink()
        yield
        logging.shutdown()

    def test_log_created_and_contains_info(self, tmp_path):
        logger = logging.getLogger()
        handler = logging.FileHandler(self.log_path)
        logger.addHandler(handler)
        logger.setLevel(logging.INFO)

        logger.info("Test log entry")

        assert self.log_path.exists(), "Log file was not created"
        log_content = self.log_path.read_text()
        assert "Test log entry" in log_content, f"Log content not found. Content: {log_content}"

        logging.shutdown()

        self.log_path.unlink()

