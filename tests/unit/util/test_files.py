from unittest import TestCase
import tempfile

from mock import patch
from rediting.util import files


class TestFiles(TestCase):

    def test_nonblank_lines(self):
        test_file = tempfile.NamedTemporaryFile(delete=False)
        with test_file:
            test_file.write("First line\n")
            test_file.write("Second line\n")
            test_file.write("\n")
            test_file.write("Third line\n")

        lines = list(files.nonblank_lines(open(test_file.name)))
        self.assertItemsEqual(lines, ["First line", "Second line", "Third line"])

    def test_build_seqdict_1_sequence(self):
        test_file = tempfile.NamedTemporaryFile(delete=False)
        with test_file:
            test_file.write(">blah\n")
            test_file.write("some value")

        input_dict = {}
        files.build_seqdict(test_file.name, input_dict)
        self.assertEqual(input_dict, {"blah": "some value"})

    @patch.object(files, "nonblank_lines")
    def test_build_seqdict_multi_sequence(self, mock_non_blank):
        mock_non_blank.return_value = [">blah", "some value", "", ">foo", "bar"]

        test_file = tempfile.NamedTemporaryFile(delete=False)
        input_dict = {}
        files.build_seqdict(test_file.name, input_dict)
        self.assertEqual(input_dict, {"blah": "some value", "foo": "bar"})
