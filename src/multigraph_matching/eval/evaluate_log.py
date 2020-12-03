#!/usr/bin/python
import re
import matplotlib.pyplot as plt
import sys
import os
import os.path
from subprocess import call
from collections import namedtuple
time_series_element = namedtuple("time_series", "value time")
time_series_concatenation_element = namedtuple("time_series", "time_series_number value time")

methods = ['tightening_mcf', 'mcf']
#methods = ['tightening_mp', 'tightening_mcf', 'mp', 'mcf']
line_colors = ['r', 'g', 'b', 'c', 'm', 'y', 'k']

def tex_escape(text):
   """
      :param text: a plain text message
       :return: the message escaped to appear correctly in LaTeX
   """
   conv = {
      '&': r'\&',
      '%': r'\%',
      '$': r'\$',
      '#': r'\#',
      '_': r'\_',
      '{': r'\{',
      '}': r'\}',
      '~': r'\textasciitilde{}',
      '^': r'\^{}',
      '\\': r'\textbackslash{}',
      '<': r'\textless{}',
      '>': r'\textgreater{}',
   }
   regex = re.compile('|'.join(re.escape(unicode(key)) for key in sorted(conv.keys(), key = lambda item: - len(item)))) 
   return regex.sub(lambda match: conv[match.group()], text)

def read_log(input_file):
   lower_bounds = []
   upper_bounds = []

   for line in open(input_file,'r'):
      if line.startswith("iteration = "):
         m = re.match("iteration = (\d+), lower bound = (-?\d+\.?\d*), time elapsed = (\d+\.?\d*)s", line)
         if m:
            iteration = m.group(1)
            lower_bound = m.group(2)
            time_elapsed = m.group(3)
            lower_bounds.append(time_series_element(lower_bound, time_elapsed))

         m = re.match("iteration = (\d+), upper bound = (-?\d+\.?\d*), time elapsed = (\d+\.?\d*)s", line)
         if m:
            iteration = m.group(1)
            upper_bound = m.group(2)
            time_elapsed = m.group(3)
            upper_bounds.append(time_series_element(upper_bound, time_elapsed))

         m = re.match("iteration = (\d+), lower bound = (-?\d+\.?\d*), upper bound = (-?\d+\.?\d*), time elapsed = (\d+\.?\d*)s", line)
         if m:
            iteration = m.group(1)
            lower_bound = m.group(2)
            upper_bound = m.group(3)
            time_elapsed = m.group(4)
            lower_bounds.append(time_series_element(lower_bound, time_elapsed))
            upper_bounds.append(time_series_element(upper_bound, time_elapsed))

   return lower_bounds, upper_bounds

def average_time_series(time_series):
   no_time_series = len(time_series)
   time_series_concatenate = []
   for i in range(0,no_time_series):
      for elem in time_series[i]:
         time_series_concatenate.append(time_series_concatenation_element (i, elem.value, elem.time))
   sorted_ts = sorted(time_series_concatenate, key=lambda x: x.time)

   initial_index = 0
   index = [0] * no_time_series
   no_covered = 0
   first_general_index = 0
   for i in range(0, len(sorted_ts)):
      time_series_no = sorted_ts[i].time_series_number
      if index[i] == 0:
         index[i] = i
         no_covered += 1
         if no_covered == no_time_series:
            first_general_index = i
            break

   current_average = 0
   for i in index:
      current_average += float(sorted_ts[i].value)

   averaged_ts = [time_series_element(current_average / float(no_time_series), sorted_ts[first_general_index].time)]
   for i in range(first_general_index, len(sorted_ts)):
      time_series_number = sorted_ts[i].time_series_number
      current_average -= float(sorted_ts[index[time_series_number]].value)
      index[time_series_number] = i
      current_average += float(sorted_ts[i].value)
      averaged_ts.append(time_series_element(current_average / float(no_time_series), sorted_ts[i].time)) 

   return averaged_ts

def average_logs(method_name, line_color, input_files):
   lower_bounds = []
   upper_bounds = []
   for input_file in input_files:
      lb, ub = read_log(input_file)
      lower_bounds.append(lb)
      upper_bounds.append(ub)

   averaged_lower_bounds = average_time_series(lower_bounds)
   averaged_upper_bounds = average_time_series(upper_bounds)

   x_axis = []
   y_axis = []
   for x in averaged_lower_bounds:
      x_axis.append(float(x.time))
      y_axis.append(float(x.value))
   plt.plot(x_axis,y_axis, label=method_name, c=line_color)
   plt.legend()

   x_axis = []
   y_axis = []
   for x in averaged_upper_bounds:
      x_axis.append(float(x.time))
      y_axis.append(float(x.value))
   plt.plot(x_axis,y_axis, c=line_color, dashes=[6, 2])
   plt.legend()

   plt.ylabel('average energy/lower bounds')
   plt.xlabel('time (s)')


def average_results(input_files):
   lower_bounds = []
   upper_bounds = []
   times = []
   get_values = lambda l: map(lambda x: float(x.value), l)
   get_times = lambda l: map(lambda x: float(x.time), l)
   for input_file in input_files:
      lb, ub = read_log(input_file)
      lower_bounds.append(max(get_values(lb)))
      upper_bounds.append(min(get_values(ub)))
      times.append(max(get_times(lb) + get_times(ub)))

   avg_lb = sum(lower_bounds) / float(len(lower_bounds))
   avg_ub = sum(upper_bounds) / float(len(upper_bounds))
   avg_time = sum(times) / float(len(times))
   return avg_lb, avg_ub, avg_time

def dataset_name(input_folder):
   if input_folder.endswith('/'): input_folder = input_folder[:-1]
   return  input_folder.split('/')[-1]

def compute_primal_dual_table(methods, input_folders):
   latex = '\\begin{table}\n\\begin{tabular}{ll' + 'c' * len(methods) +'}\n'
   methods = map(tex_escape, methods)
   latex += 'dataset && ' + ' & '.join(methods) + ' \\\\ \\hline\n'
   for input_folder in input_folders:
      latex += '\multirow{3}{*}{' + tex_escape(dataset_name(input_folder)) + '}\n'
      latex += '& UB '
      for method in methods:
         lb, ub, time = average_results(input_files[i])
         latex += ' & ' + str(lb)
      latex += '\\\\\n& LB '
      for method in methods:
         lb, ub, time = average_results(input_files[i])
         latex += ' & ' + str(ub)
      latex += '\\\\\n& time(s) '
      for method in methods:
         lb, ub, time = average_results(input_files[i])
         latex += ' & ' + str(time)
      latex += '\\\\ \\hline\n' 

   latex += '\\end{tabular}\n\\end{table}\n'
   return latex

def get_log_file_names(instance_file, algorithm):
   instance_name = os.path.splitext(os.path.basename(instance))[0]
   log_file = instance_name + "_" + algorithm + "_log.txt"
   result_file = instance_name + "_" + algorithm + "_result.txt"
   return log_file, result_file

def remove_incomplete_log(instance_file, algorithm):
   # remove log and result file if only one of them is present or either of them is empty
   log_file, result_file = get_log_file_names(instance_file, algorithm)

   if not (os.path.isfile(log_file) and os.path.isfile(result_file))
      or (os.path.isfile(log_file) and os.path.getsize(log_file) == 0)
      or (os.path.isfile(result_file) and os.path.getsize(result_file) == 0):
         os.remove(log_file)
         os.remove(result_file)

def evaluate_instance(instance_file, algorithm, executable_dir, instance_dir, output_dir, solver_options):
   full_algorithm_path = os.path.join(executable_dir, algorithm)
   instance_file = os.path.join(instance_dir, instance_file)
   output_dir = os.path.join(output_dir,  os.path.dirname(instance_file))
   instance_name = os.path.splitext(os.path.basename(instance))[0]
   log_file, result_file = get_log_file_names(instance_file, algorithm)
   log_file = os.path.join(output_dir, log_file)
   result_file = os.path.join(output_dir, result_file)

   create_dir_cmd = "mkdir -p " + output_dir 
   call(create_dir_cmd, shell=True)

  if not (os.path.isfile(log_file) or os.path.isfile(result_file)):
     print("optimize " + instance_name + " with " + algorithm)
     run_alg_cmd = "export OMP_NUM_THREADS=1; " + full_algorithm_path + " -i " + instance_file + " -o " + result_file + " " + solver_options + " 2>&1 > " + log_file
     call(run_alg_cmd, shell=True)
  else:
     print("skip optimization of " + instance + " with algorithm " + algorithm)

def run_experiments(instance_list, algorithms, executable_dir, instance_dir, output_dir, solver_options):
   for instance in instance_list:
      for algorithm in algorithms:
         evaluate_instance(instance, algorithm, executable_dir, instance_dir, output_dir, solver_options)

def create_plots(paths, methods):
   for path in paths:
      for (dirpath, dirnames, filenames) in os.walk(path):
         input_files = [[] for i in range(len(methods))]
         for f in filenames:
            for method_index in range(0, len(methods)):
               if f.endswith(methods[method_index] + "_log.txt"):
                  input_file = os.path.join(dirpath,f)
                  input_files[method_index].append(input_file)
                  break

      #plt.show()
      if path.endswith('/'): path = path[:-1]
      figure_name = path.split('/')[-1]
      plt.savefig( figure_name + '.pdf')

def create_table(parths, methods):
   table_latex = compute_primal_dual_table(methods, paths)
   table_latex = '\documentclass[preview]{standalone}\n\\usepackage{multirow}\n\\begin{document}\n' + table_latex + '\end{document}'
   open('results_table.txt', 'w').write(table_latex)
   call('pdflatex results_table.txt', shell=True)

