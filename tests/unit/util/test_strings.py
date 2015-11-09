from unittest import TestCase

from rediting.util import strings


class TestString(TestCase):

    def test_gulp_long_string(self):
        long_string = "This is a very long string"
        sub_string = strings.gulp(long_string, 0, 7)
        self.assertEqual(sub_string, "This is")

    def test_gulp_long_string_0_length(self):
        long_string = "This is a very long string"
        sub_string = strings.gulp(long_string, 0, 0)
        self.assertEqual(sub_string, "")

    def test_gulp_long_string_reverse_indices(self):
        long_string = "This is a very long string"
        sub_string = strings.gulp(long_string, 7, 0)
        self.assertEqual(sub_string, "")

    def test_gulp_index_largeer_than_string(self):
        long_string = "Some string"
        sub_string = strings.gulp(long_string, 0, 20)
        self.assertEqual(sub_string, "Some string")

    def test_gulp_nothingggg(self):
        sub_string = strings.gulp("", 0,0)
        self.assertEqual(sub_string, "")

    def test_sanitize(self):
        really_dirty_string = "This-strings-is-really-dirty"
        new_string = strings.sanitize(really_dirty_string)
        self.assertEqual(new_string, "Thisstringsisreallydirty")

    def test_sanitize_all_hyphen(self):
        all_hyphens = "---"
        new_string = strings.sanitize(all_hyphens)
        self.assertEqual(new_string, "")
