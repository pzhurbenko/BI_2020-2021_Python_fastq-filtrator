import unittest
from filter_fastq import *
import os.path

class TestParsing(unittest.TestCase):
    def setUp(self):
        self.output_base_name = 'file_name'
        self.file = '/home/pzhurb/PycharmProjects/pythonProject/DNA/fastq_example'
        self.four_lines_1 = ['@name1', 'ATATATATATATATATA', '+', 'GDEFG@DFFGGGGDEFG']
        self.four_lines_2 = ['@3\n',
                             'GGGGGGGGGGGGGGGGGGGGGGGGGGGGGCCCCGGGGGGGGGGGGGGGGGCGGGGGGGCCGGGGGGGGGGG\n',
                             '+\n',
                             'CCCFFFFFHHHHHJJJJJJJJJJFFHIJJJJJJJJJJJJJJJJJJJJJJJIJHHHHHHFDEDF;AEEEEEE\n']
        self.four_lines_3 = ['@4\n',
                             'AAAAAGCGCGCGCGCGCGAAAAAAAAAAAGCAAAAAAACAAGAAACAACAAGGGCAGAAAAAAAAAAAAAA\n',
                             '+\n',
                             'CCCFFFFFHHHHHJHIIJIIIIJJJJJJGIJJJJJIJJIIGHJJJJJIIJJDHFFFFFEDACDDDCDDDDC\n']

        self.keep_filtered = True
        self.not_keep_filtered = False
        self.gc_bounds_low = 0
        self.gc_bounds_low_ziro = 0
        self.gc_bounds_low_50 = 50
        self.gc_bounds_low_100 = 100
        self.gc_bounds_high = 100
        self.gc_bounds_high_ziro = 0
        self.gc_bounds_high_50 = 50
        self.gc_bounds_high_100 = 100
        self.min_length = 0
        self.min_length_20 = 20

        self.test_1 = ['python', 'filter_fastq.py', '--min_length', '25', '--gc_bounds', '55', '65', '--keep_filtered',
                       '--output_base_name', 'output', 'file.fastq'] # correct arguments
        self.test_2 = ['python', 'filter_fastq.py'],  # less arguments
        self.test_3 = ['python', 'filter_fastq.py', '--min_lenthg', '25', 'file.fasta']  # mistake
        self.test_4 = ['abc', 'file.fastq']  # for check if no arguments
        self.test_5 = ['--min_length', '25.5']  # incorrect/not defined value min_length
        self.test_6 = ['--gc_bounds', '55.5']  # incorrect/not defined value gc_bounds
        self.test_7 = ['--gc_bounds', '55', '45']  # incorrect order gc_bounds
        self.test_8 = ['--gc_bounds', '55', 'abc']  # one value gc_bounds
        self.test_9 = ['--output_base_name', 'file.fastq']  # output_base_name is not specified


    def test_file_name(self):
        self.assertEqual(file_name(self.test_1), 'file.fastq')
        with self.assertRaises(ValueError):
            file_name(self.test_3)

    def test_min_length(self):
        self.assertEqual((min_length(self.test_1),
                          min_length(self.test_4)), (25, 0))
        with self.assertRaises(ValueError):
            min_length(self.test_5)

    def test_keep_filtered(self):
        self.assertEqual((keep_filtered(self.test_1),
                          keep_filtered(self.test_4)), (True, False))

    def test_gc_bounds(self):
        self.assertEqual((gc_bounds(self.test_1),
                          gc_bounds(self.test_4),
                          gc_bounds(self.test_8)),
                         ((55, 65), (0, 100), (55, 100)))
        with self.assertRaises(ValueError):
            gc_bounds(self.test_6)
        with self.assertRaises(ValueError):
            gc_bounds(self.test_7)

    def test_output_base_name(self):
        self.assertEqual((output_base_name(self.test_1),
                          output_base_name(self.test_4)), ('output', 'file'))
        with self.assertRaises(ValueError):
            output_base_name(self.test_9)

    def test_output_names(self):
        self.out_failed, self.out_passed = output_names(self.output_base_name)
        self.assertEqual(self.out_failed, 'file_name__failed.fastq')
        self.assertEqual(self.out_passed, 'file_name__passed.fastq')

    def test_open_files(self):
        open_files('file_name__passed.fastq', 'file_name__failed.fastq', self.keep_filtered)
        check_file = os.path.exists('file_name__failed.fastq')  # True
        self.assertEqual(check_file, True)

        size_passed = os.path.getsize('file_name__passed.fastq') > 0
        self.assertEqual(size_passed, False)

        size_failed = os.path.getsize('file_name__failed.fastq') > 0
        self.assertEqual(size_failed, False)

        os.remove('file_name__passed.fastq')
        os.remove('file_name__failed.fastq')

        open_files('file_name__passed.fastq', 'file_name__failed.fastq', self.not_keep_filtered)
        check_file = os.path.exists('file_name__failed.fastq')  # True
        self.assertEqual(check_file, False)

        os.remove('file_name__passed.fastq')

    def test_write_file(self):
        open_files('file_name__passed.fastq', 'file_name__failed.fastq', self.keep_filtered)
        filter(self.four_lines_1, self.keep_filtered, 'file_name__failed.fastq')
        size_failed = os.path.getsize('file_name__failed.fastq') > 0
        self.assertEqual(size_failed, True)

        write_file('file_name__passed.fastq', self.four_lines_1)
        size_passed = os.path.getsize('file_name__failed.fastq') > 0
        self.assertEqual(size_passed, True)

        os.remove('file_name__passed.fastq')
        os.remove('file_name__failed.fastq')

    def test_fastq_parse(self):
        four_lines_iterator = fastq_parse(self.file)
        four_lines = fastq_parse_reader(four_lines_iterator)
        self.assertEqual(four_lines, self.four_lines_1)

    def test_gc_count(self):
        gc_2 = gc_count(self.four_lines_2[1][:-1])
        self.assertEqual(gc_2, 100)
        gc_3 = gc_count(self.four_lines_3[1][:-1])
        self.assertEqual(gc_3, 33)
        gc_1 = gc_count(self.four_lines_1[1][:-1])
        self.assertEqual(gc_1, 0)

    def test_filter_gc_content(self):
        gc_2_yes = filter_gc_content(self.four_lines_2, self.gc_bounds_low_ziro, self.gc_bounds_high)
        self.assertEqual(gc_2_yes, True)
        gc_2_not = filter_gc_content(self.four_lines_2, self.gc_bounds_low_ziro, self.gc_bounds_high_50)
        self.assertEqual(gc_2_not, False)
        gc_3_not = filter_gc_content(self.four_lines_3, self.gc_bounds_low_50, self.gc_bounds_high)
        self.assertEqual(gc_3_not, False)
        gc_3_yes = filter_gc_content(self.four_lines_3, self.gc_bounds_low_ziro, self.gc_bounds_high_ziro)
        self.assertEqual(gc_3_yes, False)


if __name__ == "__main__":
    unittest.main()
