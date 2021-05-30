import sys

inp = sys.argv[1:]

def min_length(inp):
    if '--min_length' in inp:
        length_index = inp.index('--min_length')
        try:
            length_value = int(inp[length_index + 1])
        except ValueError:
            raise ValueError('Вы не указали числовое значение min_length')
    else:
        length_value = 0
    return length_value

def keep_filtered(inp):
    if '--keep_filtered' in inp:
        return True
    else:
        return False

def gc_bounds(inp):
    if '--gc_bounds' in inp:
        gc_index = inp.index('--gc_bounds')
        try:
            gc_value_min = int(inp[gc_index + 1])
        except ValueError:
            raise ValueError('Вы не указали числовое значение gc_bounds')
        try:
            gc_value_max = int(inp[gc_index + 2])
        except ValueError:
            gc_value_max = 100
        if gc_value_min > gc_value_max:
            raise ValueError('Укажите сначала нижнюю границу, затем верхнюю')
    else:
        gc_value_min = 0
        gc_value_max = 100
    return gc_value_min, gc_value_max

def gc_count(seq):
    if len(seq) > 0:
        GC = int((seq.count('G') + seq.count('C')) / len(seq) * 100)
        return GC
    else:
        return None

def fastq_parse(input_file):
    with open(input_file, 'r') as fq_in:
        lines = []
        for line in fq_in:
            lines.append(line.rstrip())
            if len(lines) == 4:
                yield lines
                lines = []

def fastq_parse_reader(lines):
    try:
        four_lines = next(lines)
        return four_lines
    except StopIteration:
        return None

def output_base_name(inp):
    if '--output_base_name' in inp:
        base_index = inp.index('--output_base_name')
        base_name = inp[base_index + 1]
        if inp[base_index + 1].endswith('.fastq') or inp[base_index + 1].startswith('--'):
            raise ValueError('Output base name is set, but not specified!')
    else:
        base_name = inp[-1].split('/')[-1].replace('.fastq', '')
    return base_name

def output_names(base_name):
    out_passed = str(base_name + '__passed.fastq')
    out_failed = str(base_name + '__failed.fastq')
    return (out_failed, out_passed)

def file_name(inp):
    if inp[-1].endswith('.fastq'):
        return inp[-1]
    else:
        raise ValueError('Последний аргумент должен являться fastq файлом')

def open_files(out_passed, out_failed, keep_filtered):
    with open(out_passed, 'w') as f:
        pass
    if keep_filtered == True:
        with open(out_failed, 'w') as f:
            pass

def write_file(out_passed, four_lines):
    with open(out_passed, 'a') as fq_in:
        for i in four_lines:
            fq_in.writelines("%s\n" % i)

def filter(four_lines, keep_filtered, out_failed):
    if keep_filtered == True:
        with open(out_failed, 'a') as fq_out:
            for i in four_lines:
                fq_out.writelines("%s\n" % i)

def filter_gc_content(four_lines, gc_value_min, gc_value_max):
    gc_content = gc_count(four_lines[1])
    if gc_content < gc_value_min or gc_content > gc_value_max:
        return True
    return False

if __name__ == '__main__':
    file = file_name(inp)
    min_length = min_length(inp)
    keep_filtered = keep_filtered(inp)
    gc_value_min, gc_value_max = gc_bounds(inp)
    output_base_name = output_base_name(inp)
    out_failed, out_passed = output_names(output_base_name)
    open_files(out_passed, out_failed, keep_filtered)
    four_lines_iterator = fastq_parse(file)
    while True:
        four_lines = fastq_parse_reader(four_lines_iterator)
        if four_lines is None:
            break
        if len(four_lines[1]) < min_length:
            if keep_filtered == True:
                filter(four_lines, keep_filtered, out_failed)
            continue
        if gc_value_min is not None:
            if filter_gc_content(four_lines, gc_value_min, gc_value_max) == False:
                write_file(out_passed, four_lines)
            else:
                if keep_filtered == True:
                    filter(four_lines, keep_filtered, out_failed)
                continue
        else:
            filter(four_lines, keep_filtered, out_failed)
            write_file(out_passed, four_lines)
            continue
    print('Успешно')
