import sys

inp = list(sys.argv[1:])

# Определяем параметры:
input_file = sys.argv[-1]
#
if '--min_length' in inp:
    length_index = inp.index('--min_length')
    length_value = int(inp[length_index + 1])
else:
    length_value = 0
#
if '--gc_bounds' in inp:
    gc_index = inp.index('--gc_bounds')
    gc_value_min = int(inp[gc_index + 1])
    gc_value_max = 100
    if inp[gc_index + 2].isnumeric():
        gc_value_max = int(inp[gc_index + 2])
else:
    gc_value_min = 0
    gc_value_max = 100
#
if '--keep_filtered' in inp:
    filtered = True
else:
    filtered = False
#
if '--output_base_name' in inp:
    base_index = inp.index('--output_base_name')
    base_name = inp[base_index + 1]
else:
    base_name = input_file[0: -6]

# Функция для подсчета GC состава:
def gc_count(seq):
    global gc_value_min
    global gc_value_max
    GC = int((seq.count('G') + seq.count('C')) / len(seq) * 100)
    return GC

# Основная функция для обработки fastq:
def fastq_filter(input_file, base_name, length_value, filtered, gc_value_min, gc_value_max):

    # Если нужно записать не прошедшии фильтрацию риды в failed:
    if filtered == True:
        with open(input_file, 'r') as fq_in, open('output_' + base_name + '__passed.fastq', 'a') as fq_out, open('output_' + base_name + '__failed.fastq', 'a') as fq_failed:
            # Набираем по 4 строки в list
            lines = []
            count = 0
            for line in fq_in:
                lines.append(line.rstrip())
                if len(lines) == 4:
                    count += 1
                    seq = lines[1]

                    # Фильтрация по GC составу
                    GC = gc_count(seq)
                    if  GC >= gc_value_min and GC <= gc_value_max:

                         # Фильтрация по длине
                         if len(seq) >= length_value:
                             lines[0] = '@' + base_name + '_' + str(count)
                             fq_out.writelines("%s\n" % i for i in lines)
                             lines = []

                         # Запись в failed
                         else:
                             lines[0] = '@' + base_name + '_' + str(count)
                             fq_failed.writelines("%s\n" % i for i in lines)
                             lines = []

                    # Запись в failed
                    else:
                        lines[0] = '@' + base_name + '_' + str(count)
                        fq_failed.writelines("%s\n" % i for i in lines)
                        lines = []

    # Если НЕ нужно записать не прошедшии фильтрацию риды в failed:  (как сделать оптимальнее?!)
    else:
        with open(input_file, 'r') as fq_in, open('output_' + base_name + '__passed.fastq', 'a') as fq_out:
            lines = []
            count = 0
            for line in fq_in:
                lines.append(line.rstrip())
                if len(lines) == 4:
                     count += 1
                     seq = lines[1]

                     GC = gc_count(seq)
                     if  GC >= gc_value_min and GC <= gc_value_max:
                          if len(seq) >= length_value:
                              lines[0] = '@' + base_name + '_' + str(count)
                              fq_out.writelines("%s\n" % i for i in lines)
                     lines = []

# запускаем функцию
fastq_filter(input_file, base_name, length_value, filtered, gc_value_min, gc_value_max)


