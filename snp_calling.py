import subprocess
import time
import os
import csv

# Пути для программ
FASTQC = "fastqc * -t 72" 
trimmomatic_jar = '/home/ann/Tools/Trimmomatic-main/trimmomatic-0.39.jar'
GATK = '/home/ann/Tools/gatk-4.1.8.1/gatk-package-4.1.8.1-local.jar'
FreeBayes = 'freebayes'
picard = 'java -jar /home/ann/Tools/picard.jar'


# Путь к CSV-файлу для логирования времени выполнения команд
LOG_FILE = '/home/ann/Variant_calling/bcf_bwamem2.csv'

# Создаем папку, если её нет
os.makedirs(os.path.dirname(LOG_FILE), exist_ok=True)

# Функция для выполнения команды и записи времени выполнения в CSV-файл
def run_command1(command):
    start_time = time.time()  # Запоминаем время начала
    print(f"Running command: {command}")
    try:
        subprocess.run(command, shell=True, check=True)
    except subprocess.CalledProcessError as e:
        print(f"Error occurred: {e}")
    elapsed_time = time.time() - start_time  # Вычисляем время выполнения
    print(f"Command finished in {elapsed_time:.2f} seconds")
    
    # Запись данных в CSV-файл
    file_exists = os.path.isfile(LOG_FILE)  # Проверяем, существует ли файл

    with open(LOG_FILE, 'a', newline='') as log_file:
        writer = csv.writer(log_file)
        # Добавляем заголовки, если файл создается впервые
        if not file_exists:
            writer.writerow(["Command", "Elapsed Time (seconds)"])
        # Записываем команду и время выполнения
        writer.writerow([command, f"{elapsed_time:.2f}"])
    return elapsed_time

def run_command(command, folder_name):
    start_time = time.time()  # Запоминаем время начала
    print(f"Running command for {folder_name}: {command}")
    
    try:
        subprocess.run(command, shell=True, check=True)
    except subprocess.CalledProcessError as e:
        print(f"Error occurred: {e}")
    
    elapsed_time = time.time() - start_time  # Вычисляем время выполнения
    print(f"Command for {folder_name} finished in {elapsed_time:.2f} seconds")
    
    # Запись данных в CSV-файл
    file_exists = os.path.isfile(LOG_FILE)  # Проверяем, существует ли файл
    with open(LOG_FILE, 'a', newline='') as log_file:
        writer = csv.writer(log_file)
        if not file_exists:
            writer.writerow(["Folder Name", "Command", "Elapsed Time (seconds)"])
        writer.writerow([folder_name, command, f"{elapsed_time:.2f}"])
    
    return elapsed_time

# Функция для контроля качества (FastQC) файлов .fastq
def quality_control(fastq_files_path, output_dir):
    original_dir = os.getcwd()  # Сохраняем текущую рабочую директорию
    command = f"cd {fastq_files_path} && {FASTQC} -o {output_dir} && cd -"  # Добавляем cd - для возврата в исходную директорию
    run_command1(command)
    os.chdir(original_dir)  # Возвращаемся в исходную директорию после выполнения команды
   #print(f"Current directory after cd: {os.getcwd()}")  # для проверки

# Фильтрация
def trim_reads(fastq_files_path, output_dir, trimmomatic_jar):
    total_elapsed_time = 0
    processed_pairs = set()  # Множество для хранения обработанных пар файлов

    # Перебор всех файлов в указанной директории
    for file_name in os.listdir(fastq_files_path):
        if file_name.endswith("_R1_001.fastq.gz"):
            folder_name = file_name.split("_R1")[0]  # Извлечение части имени файла перед _R1

            input_file_R1 = os.path.join(fastq_files_path, file_name)  # Полный путь к файлу R1
            input_file_R2 = os.path.join(fastq_files_path, file_name.replace("_R1_001.fastq.gz", "_R2_001.fastq.gz"))
            file_pair = (input_file_R1, input_file_R2)
            if file_pair not in processed_pairs:  # Проверка, была ли пара файлов уже обработана
                output_file_R1_paired = os.path.join(output_dir, f"{folder_name}_paired_R1.fastq.gz")
                output_file_R1_unpaired = os.path.join(output_dir, f"{folder_name}_unpaired_R1.fastq.gz")
                output_file_R2_paired = os.path.join(output_dir, f"{folder_name}_paired_R2.fastq.gz")
                output_file_R2_unpaired = os.path.join(output_dir, f"{folder_name}_unpaired_R2.fastq.gz")
                
                # Команда для запуска Trimmomatic
                command = f"{trimmomatic_jar} PE  {input_file_R1} {input_file_R2} {output_file_R1_paired} {output_file_R1_unpaired} {output_file_R2_paired} {output_file_R2_unpaired} -threads 72 ILLUMINACLIP:/home/ann/Tools/Trimmomatic-main/adapters/TruSeq3-PE.fa:2:30:10:2:True"
                elapsed_time = run_command1(command)
                total_elapsed_time += elapsed_time  # Суммирование общего времени
                processed_pairs.add(file_pair)  # Добавление пары файлов в множество processed_pairs  

    print(f"Total time for processing: {total_elapsed_time:.2f} seconds")  # Вывод общего времени обработки 

    
# Выравнивание

# bwa mem
def bwa_mem(input_files, output_dir, index_ref_file):
    total_elapsed_time = 0
    processed_pairs = set()  # Множество для хранения обработанных пар файлов

    for file_name in os.listdir(input_files):
        if file_name.endswith("_paired_R1.fastq.gz"):
            folder_name = file_name.split("_paired")[0]  # Извлечение части имени файла перед _paired
            # Создание пути для файла .sam 
            output_file_sam = os.path.join(output_dir, f"{folder_name}.sam")
            os.makedirs(output_dir, exist_ok=True)

            input_file_R1 = os.path.join(input_files, file_name) #Создание полного пути к файлу R1
            input_file_R2 = os.path.join(input_files, file_name.replace("_R1.fastq.gz", "_R2.fastq.gz"))
            file_pair = (input_file_R1, input_file_R2)

            if file_pair not in processed_pairs: #Проверка, была ли эта пара файлов уже обработана
                command = (
                    f"bwa mem -t 72 {index_ref_file} {input_file_R1} {input_file_R2} > {output_file_sam} "
                )
                elapsed_time = run_command(command, folder_name)
                total_elapsed_time += elapsed_time
                processed_pairs.add(file_pair) #Добавление пары файлов в множество processed_pairs

    print(f"Total time for processing: {total_elapsed_time:.2f} seconds") 
    return total_elapsed_time
    
def bwa_mem2(input_files, output_dir, index_ref_file):
    total_elapsed_time = 0
    processed_pairs = set()  # Множество для хранения обработанных пар файлов

    for file_name in os.listdir(input_files):
        if file_name.endswith("_paired_R1.fastq.gz"):
            folder_name = file_name.split("_paired")[0]  # Извлечение части имени файла перед _paired
            # Создание пути для файла .sam 
            output_file_sam = os.path.join(output_dir, f"{folder_name}.sam")
            os.makedirs(output_dir, exist_ok=True)

            input_file_R1 = os.path.join(input_files, file_name) #Создание полного пути к файлу R1
            input_file_R2 = os.path.join(input_files, file_name.replace("_R1.fastq.gz", "_R2.fastq.gz"))
            file_pair = (input_file_R1, input_file_R2)

            if file_pair not in processed_pairs: #Проверка, была ли эта пара файлов уже обработана
                command = (
                    f"bwa-mem2.sse41 mem -t 72 {index_ref_file} {input_file_R1} {input_file_R2} > {output_file_sam} "
                )

                elapsed_time = run_command(command, folder_name)
                total_elapsed_time += elapsed_time
                processed_pairs.add(file_pair) #Добавление пары файлов в множество processed_pairs

    print(f"Total time for processing: {total_elapsed_time:.2f} seconds") 
    return total_elapsed_time
    
    

# minimap2
def minimap2(input_files, output_dir, index_ref_file):
    total_elapsed_time = 0
    processed_pairs = set()  # Множество для хранения обработанных пар файлов

    for file_name in os.listdir(input_files):
        if file_name.endswith("_paired_R1.fastq.gz"):
            folder_name = file_name.split("_paired")[0]  # Извлечение части имени файла перед _paired
            # Создание пути для файла .sam 
            output_file_sam = os.path.join(output_dir, f"{folder_name}.sam")
            os.makedirs(output_dir, exist_ok=True)

            input_file_R1 = os.path.join(input_files, file_name) #Создание полного пути к файлу R1
            input_file_R2 = os.path.join(input_files, file_name.replace("_R1.fastq.gz", "_R2.fastq.gz"))
            file_pair = (input_file_R1, input_file_R2)

            if file_pair not in processed_pairs: #Проверка, была ли эта пара файлов уже обработана
                command = (
                    f"minimap2 -t 72 -a {index_ref_file} {input_file_R1} {input_file_R2} > {output_file_sam} "
                )

                elapsed_time = run_command(command, folder_name)
                total_elapsed_time += elapsed_time
                processed_pairs.add(file_pair) #Добавление пары файлов в множество processed_pairs

    print(f"Total time for processing: {total_elapsed_time:.2f} seconds") 
    return total_elapsed_time
    
# bowtie2
def bowtie2(input_files, output_dir, index_ref_file):
    total_elapsed_time = 0
    processed_pairs = set()  # Множество для хранения обработанных пар файлов

    for file_name in os.listdir(input_files):
        if file_name.endswith("_paired_R1.fastq.gz"):
            folder_name = file_name.split("_paired")[0]  # Извлечение части имени файла перед _paired
            # Создание пути для файла .sam 
            output_file_sam = os.path.join(output_dir, f"{folder_name}.sam")
            os.makedirs(output_dir, exist_ok=True)

            input_file_R1 = os.path.join(input_files, file_name) #Создание полного пути к файлу R1
            input_file_R2 = os.path.join(input_files, file_name.replace("_R1.fastq.gz", "_R2.fastq.gz"))
            file_pair = (input_file_R1, input_file_R2)

            if file_pair not in processed_pairs: #Проверка, была ли эта пара файлов уже обработана
                command = (
                    f"bowtie2 -p 72 -x {index_ref_file} -1 {input_file_R1} -2 {input_file_R2} -S {output_file_sam} "
                )

                elapsed_time = run_command(command, folder_name)
                total_elapsed_time += elapsed_time
                processed_pairs.add(file_pair) #Добавление пары файлов в множество processed_pairs

    print(f"Total time for processing: {total_elapsed_time:.2f} seconds") 
    
    return total_elapsed_time

#sam->bam
def sam_to_bam(sam_dir):

    for file_name in os.listdir(sam_dir):
        if file_name.endswith(".sam"):
            # Полные пути к входному SAM и выходному BAM
            input_file_sam = os.path.join(sam_dir, file_name)
            output_file_bam = os.path.join(sam_dir, file_name.replace(".sam", ".bam"))
            sorted_bam = os.path.join(sam_dir, file_name.replace(".sam", "_sorted.bam"))
            
            # Команда для преобразования SAM в BAM
            command = ( 
            f"samtools view -@ 72 -Sb {input_file_sam} > {output_file_bam} &&"
            f"samtools sort -@ 72 {output_file_bam} -o {sorted_bam}")
            
            # Выполняем команду
            run_command1(command)


def run_quality(command):
    # Выполняем команду и захватываем её вывод
    result = subprocess.run(command, shell=True, capture_output=True, text=True)
    return result.stdout.strip()  # Возвращаем вывод без лишних пробелов и символов новой строки
               
      
def alignment_quality(bam_dir):
    output_file = os.path.join(bam_dir, "alignment_quality_results.csv")
    
    # Открываем файл для записи в формате CSV с разделителем ";"
    with open(output_file, 'w', newline='') as out_f:
        writer = csv.writer(out_f, delimiter=';')  # Используем точку с запятой как разделитель
        
        # Записываем заголовок таблицы (названия параметров)
        writer.writerow(["Файл", "Процент покрытия", "Средняя глубина", "Количество выравненных прочтений", "Среднее MQ", "Медиана MQ"])
        
        # Перебираем все BAM файлы в указанной директории
        for file_name in os.listdir(bam_dir):
            if file_name.endswith("_sorted.bam"):
                # Полный путь к входному BAM файлу
                input_file = os.path.join(bam_dir, file_name)

                # Команды для выполнения
                command_depth = f"export LC_ALL=C && samtools depth -a -b '/путь/к/Exome.bed' {input_file} | awk '{{if($3 >= 8) covered++}} END {{print (covered/NR)*100}}'"
                command_avg_depth = f"export LC_ALL=C && samtools depth {input_file} | awk '{{sum+=$3}} END {{print sum/NR}}'"
                command_flagstat = f"export LC_ALL=C && samtools flagstat {input_file} | grep 'mapped (' | head -n 1 | awk '{{print $1, $5\")\"}}'"
                command_quality_avg = f"export LC_ALL=C && samtools view {input_file} | cut -f 5 | awk '{{sum += $1; count++}} END {{print sum / count}}'"
                command_quality_median = f"export LC_ALL=C && samtools view {input_file} | cut -f 5 | sort -n | awk '{{values[NR] = $1}} END {{if (NR % 2) {{print values[(NR + 1) / 2]}} else {{print (values[NR / 2] + values[(NR / 2) + 1]) / 2}}}}'"


                # Выполнение команд и получение результатов
                result_depth = run_quality(command_depth).strip()
                result_avg_depth = run_quality(command_avg_depth).strip()
                result_flagstat = run_quality(command_flagstat).strip()
                result_quality_avg = run_quality(command_quality_avg).strip()
                result_quality_median = run_quality(command_quality_median).strip()

                # Записываем результаты в формате CSV (каждое значение разделяется точкой с запятой)
                writer.writerow([file_name, result_depth, result_avg_depth, result_flagstat, result_quality_avg, result_quality_median])



#Вызов вариантов

#gatk

def gatk(input_files, output_dir, index_ref_file):
    total_elapsed_time = 0
    counter = 1

    for file_name in os.listdir(input_files):
        if file_name.endswith("_sorted.bam"):
            # Создание пути для файлов
            folder_name = file_name.split("_sorted.bam")[0]  # Извлечение части имени файла перед _sorted.bam
            input_file_bam = os.path.join(input_files, file_name) #Создание полного пути к файлу R1
            RGr_bam = os.path.join(output_dir, f"{folder_name}_RGr.bam")
            MrkDup_bam = os.path.join(output_dir, f"{folder_name}_MrkDup.bam")
            Srt_bam = os.path.join(output_dir, f"{folder_name}_Srt.bam")
            #Rcl_table = os.path.join(output_dir, f"{folder_name}_Rcl.table")
            #Rcl_bam = os.path.join(output_dir, f"{folder_name}_Rcl.bam")
            
            vcf_dir = os.path.join(output_dir, "VCF")
            os.makedirs(vcf_dir, exist_ok=True)
            final_vcf = os.path.join(vcf_dir, f"{folder_name}_final.vcf")

            os.makedirs(output_dir, exist_ok=True)

            command = (
                f" java -jar {GATK} AddOrReplaceReadGroups -I {input_file_bam}   -O {RGr_bam}   -SORT_ORDER coordinate -RGPL ILLUMINA -RGLB lib{counter}  -RGPU unit{counter}  -RGSM {folder_name} && "
                f"samtools index {RGr_bam} &&"
                f" java -jar {GATK} MarkDuplicates -I {RGr_bam}   -O {MrkDup_bam}  -M {output_dir}/{folder_name}_MrkDup_metrics.txt &&"
                f"samtools sort {MrkDup_bam} -o {Srt_bam} &&"
                f"samtools index {Srt_bam} &&"
               # f"java -jar {GATK} BaseRecalibrator -I {Srt_bam} -R {index_ref_file} --known-sites {snp_file} -O {Rcl_table} &&"
               # f"java -jar {GATK} ApplyBQSR -I {Srt_bam} -R {index_ref_file} --bqsr-recal-file {Rcl_table} -O {Rcl_bam} &&"
                f"java -Xms8g -Xmx10g -jar {GATK} HaplotypeCaller  -R {index_ref_file} -I {Srt_bam} -O {final_vcf} --native-pair-hmm-threads 72 " 
            )

            elapsed_time = run_command(command, folder_name)
            total_elapsed_time += elapsed_time
            
            counter += 1

    print(f"Total time for processing: {total_elapsed_time:.2f} seconds") 
    return total_elapsed_time


def freebayes(input_files, output_dir, index_ref_file):
    total_elapsed_time = 0
    counter = 1

    for file_name in os.listdir(input_files):
        if file_name.endswith("_sorted.bam"):
            # Создание пути для файлов
            folder_name = file_name.split("_sorted.bam")[0]  # Извлечение части имени файла перед _sorted.bam
            input_file_bam = os.path.join(input_files, file_name) #Создание полного пути к файлу R1
            RGr_bam = os.path.join(output_dir, f"{folder_name}_RGr.bam")
            MrkDup_bam = os.path.join(output_dir, f"{folder_name}_MrkDup.bam")
            Srt_bam = os.path.join(output_dir, f"{folder_name}_Srt.bam")
            
            vcf_dir = os.path.join(output_dir, "VCF")
            os.makedirs(vcf_dir, exist_ok=True)
            final_vcf = os.path.join(vcf_dir, f"{folder_name}_final.vcf")

            os.makedirs(output_dir, exist_ok=True)

            command = (
                f" {picard} AddOrReplaceReadGroups -I {input_file_bam}   -O {RGr_bam}   -SORT_ORDER coordinate -RGPL ILLUMINA -RGLB lib{counter}  -RGPU unit{counter}  -RGSM {folder_name} && "
                f"samtools index {RGr_bam} &&"
                f" {picard} MarkDuplicates -I {RGr_bam}   -O {MrkDup_bam}  -M {output_dir}/{folder_name}_MrkDup_metrics.txt &&"
                f"samtools sort {MrkDup_bam} -o {Srt_bam} &&"
                f"samtools index {Srt_bam} &&"
                f"{FreeBayes} -f {index_ref_file} {Srt_bam} > {final_vcf}"
            )

            elapsed_time = run_command(command, folder_name)
            total_elapsed_time += elapsed_time
            
            counter += 1

    print(f"Total time for processing: {total_elapsed_time:.2f} seconds") 
    return total_elapsed_time
    
def sam_bcf(input_files, output_dir, index_ref_file):
    total_elapsed_time = 0
    counter = 1

    for file_name in os.listdir(input_files):
        if file_name.endswith("_sorted.bam"):
            # Создание пути для файлов
            folder_name = file_name.split("_sorted.bam")[0]  # Извлечение части имени файла перед _sorted.bam
            input_file_bam = os.path.join(input_files, file_name) #Создание полного пути к файлу R1
            RGr_bam = os.path.join(output_dir, f"{folder_name}_RGr.bam")
            RGr_sort = os.path.join(output_dir, f"{folder_name}_sort_RGr.bam")
            MrkDup_bam = os.path.join(output_dir, f"{folder_name}_MrkDup.bam")
            Srt_bam = os.path.join(output_dir, f"{folder_name}_Srt.bam")
            file_bcf = os.path.join(output_dir, f"{folder_name}.bcf")
            
            vcf_dir = os.path.join(output_dir, "VCF")
            os.makedirs(vcf_dir, exist_ok=True)
            final_vcf = os.path.join(vcf_dir, f"{folder_name}_final.vcf")

            os.makedirs(output_dir, exist_ok=True)

            command = (
                f"samtools addreplacerg -r 'ID:group{counter}' -r 'SM:sample{counter}' -o {RGr_bam} {input_file_bam} && "
                f"samtools sort {RGr_bam} -o {RGr_sort} && "
                f"samtools index {RGr_sort} && "
                f"{picard} MarkDuplicates -I {RGr_sort}   -O {MrkDup_bam}  -M {output_dir}/{folder_name}_MrkDup_metrics.txt &&"
                f"samtools sort {MrkDup_bam} -o {Srt_bam} && "
                f"samtools index {Srt_bam} && "
                f"bcftools mpileup --threads 72 -Ou -f {index_ref_file} {Srt_bam} > {file_bcf} && "
                f"bcftools call --threads 72 -mv -o {final_vcf} {file_bcf} "
            )


            elapsed_time = run_command(command, folder_name)
            total_elapsed_time += elapsed_time
            
            counter += 1

    print(f"Total time for processing: {total_elapsed_time:.2f} seconds") 
    return total_elapsed_time

def variant_quality(vcf_dir):
    vcf_folder = os.path.join(vcf_dir, "VCF")  # Указываем папку VCF внутри директории
    output_file = os.path.join(vcf_folder, "variant_quality_results.csv")
    
    # Проверяем, существует ли директория, и если нет, создаем её
    os.makedirs(vcf_folder, exist_ok=True)

    # Открываем файл для записи в формате CSV с разделителем ";"
    with open(output_file, 'w', newline='') as out_f:
        writer = csv.writer(out_f, delimiter=';')  # Используем точку с запятой как разделитель
        
        # Записываем заголовок таблицы (названия параметров)
        writer.writerow(["Файл", "Кол-во вариантов", "Среднее качество вариантов", "Медиана качества вариантов"])
        
        # Перебираем все VCF файлы в указанной директории
        for file_name in os.listdir(vcf_folder):
            if file_name.endswith("_final.vcf"):
                # Полный путь к входному VCF файлу
                input_file = os.path.join(vcf_folder, file_name)

                # Команды для выполнения
                command_variant_count = f"export LC_ALL=C && bcftools view {input_file} | grep -v '^#' | wc -l"
                command_quality_avg = f"export LC_ALL=C && bcftools view {input_file} | awk '{{if ($1 !~ /^#/ && $6 != \".\") print $6}}' | awk '{{sum += $1; count++}} END {{print sum / count}}'"
                command_quality_median = f"export LC_ALL=C && bcftools view {input_file} | awk '{{if ($1 !~ /^#/ && $6 != \".\") print $6}}' | sort -n | awk '{{values[NR] = $1}} END {{if (NR % 2) {{print values[(NR + 1) / 2]}} else {{print (values[NR / 2] + values[(NR / 2) + 1]) / 2}}}}'"



                # Выполнение команд и получение результатов
                result_variant_count = run_quality(command_variant_count)
                result_avg_quality = run_quality(command_quality_avg)
                result_median_quality = run_quality(command_quality_median)
                
                # Записываем результаты в формате CSV (каждое значение разделяется точкой с запятой)
                writer.writerow([file_name, result_variant_count, result_avg_quality, result_median_quality])



# Основная функция
def main():
    # Пути к входным и выходным данным
    fastq_files_path = '/home/ann/Variant_calling/data/Merged_Fastq/'
    output_dir = '/home/ann/Variant_calling/fastqc2'
    trimmed_files = '/home/ann/Variant_calling/Trim'
    
    index_ref_file = '/home/ann/Variant_calling/Human_ref/GRCh37_only_chr_lin.fna'
    index_ref_bowtie = '/home/ann/Variant_calling/Human_ref/GRCh37_only_chr_lin_index'
    index_ref_bwamem2 = '/home/ann/Variant_calling/Human_ref/bwa_mem2/GRCh37_only_chr_lin.fna'
    
    out_bwa = '/home/ann/Variant_calling/bwa_mem'
    out_minimap = '/home/ann/Variant_calling/minimap2'
    out_bowtie = '/home/ann/Variant_calling/bowtie2'
    out_bwamem2 = '/home/ann/Variant_calling/bwamem2'
    
    out_gatk_bwa = '/home/ann/Variant_calling/gatk/bwa_mem'
    out_gatk_minimap2 = '/home/ann/Variant_calling/gatk/minimap2'
    out_gatk_bowtie2 = '/home/ann/Variant_calling/gatk/bowtie2'
    out_gatk_bwamem2 = '/home/ann/Variant_calling/gatk/bwamem2'
    
    out_freebayes_bwa = '/home/ann/Variant_calling/freebayes/bwa_mem'
    out_freebayes_minimap2 = '/home/ann/Variant_calling/freebayes/minimap2'
    out_freebayes_bowtie2 = '/home/ann/Variant_calling/freebayes/bowtie2'
    out_freebayes_bwamem2 = '/home/ann/Variant_calling/freebayes/bwamem2'

    out_sambcf_bwa = '/home/ann/Variant_calling/sam_and_bcf/bwa_mem'
    out_sambcf_minimap2 = '/home/ann/Variant_calling/sam_and_bcf/minimap2'
    out_sambcf_bowtie2 = '/home/ann/Variant_calling/sam_and_bcf/bowtie2'
    out_sambcf_bwamem2 = '/home/ann/Variant_calling/sam_and_bcf/bwamem2'
    

    # Запуск функций
    # Оценка качества
    #quality_control(fastq_files_path, output_dir)
    
    # Фильтрация
    #trim_reads(fastq_files_path, trimmed_files, trimmomatic_jar)
    
    #bwa_mem(fastq_files_path, out_bwa, index_ref_file)
    #minimap2(fastq_files_path, out_minimap, index_ref_file)
    #bowtie2(fastq_files_path, out_bowtie, index_ref_bowtie)
    #bwa_mem2 (fastq_files_path, out_bwamem2, index_ref_bwamem2)
        
    #sam_to_bam(out_bwa)
    #sam_to_bam(out_minimap)
    #sam_to_bam(out_bowtie)
    #sam_to_bam(out_bwamem2)
    #alignment_quality(out_bwa)
    #alignment_quality(out_minimap)
    #alignment_quality(out_bowtie)
    #alignment_quality(out_bwamem2)
    
    #gatk(out_bwa, out_gatk_bwa, index_ref_file)
    #gatk(out_minimap, out_gatk_minimap2, index_ref_file)
    #gatk(out_bowtie, out_gatk_bowtie2, index_ref_file)
    #gatk(out_bwamem2, out_gatk_bwamem2, index_ref_file)
    
    #freebayes(out_bwa, out_freebayes_bwa, index_ref_file)
    #freebayes(out_minimap, out_freebayes_minimap2, index_ref_file)
    #freebayes(out_bowtie, out_freebayes_bowtie2, index_ref_file)
    #freebayes(out_bwamem2, out_freebayes_bwamem2, index_ref_file)
    
    #sam_bcf(out_bwa, out_sambcf_bwa, index_ref_file)
    #sam_bcf(out_minimap, out_sambcf_minimap2, index_ref_file)
    #sam_bcf(out_bowtie, out_sambcf_bowtie2, index_ref_file)
    #sam_bcf(out_bwamem2, out_sambcf_bwamem2, index_ref_file)
    
    #variant_quality(out_gatk_bwa)
    #variant_quality(out_gatk_minimap2)
    #variant_quality(out_gatk_bowtie2)
    #variant_quality(out_gatk_bwamem2)
    #variant_quality(out_freebayes_bwa)
    #variant_quality(out_freebayes_minimap2)
    #variant_quality(out_freebayes_bowtie2)
    #variant_quality(out_freebayes_bwamem2)
    #variant_quality(out_sambcf_bwa)
    #variant_quality(out_sambcf_minimap2)
    #variant_quality(out_sambcf_bowtie2)
    #variant_quality(out_sambcf_bwamem2)

if __name__ == "__main__":
    main()


