#!/usr/bin/python
import re
import matplotlib
matplotlib.use('Agg') # to avoid problems with no display
import matplotlib.pyplot as plt
import sys
import os
import os.path
from subprocess import call
from collections import namedtuple

# TODO: make class out of this

time_series_element = namedtuple("time_series", "value time")
time_series_concatenation_element = namedtuple("time_series", "time_series_number value time")

line_colors = ['r', 'g', 'b', 'c', 'm', 'y', 'k']

def get_log_file_names(instance_file, algorithm, output_dir, opts_name):
   instance_dir = os.path.dirname(instance_file)
   instance_name = os.path.splitext(os.path.basename(instance_file))[0]
   ending = os.path.splitext(os.path.basename(instance_file))[1]
   if(ending not in ('.lp', '.txt', '.uai', '.hdf5', '.h5')):
       instance_name = instance_name + ending
   log_file = instance_name + "_" + algorithm + "_" + opts_name + "_log.txt"
   result_file = instance_name + "_" + algorithm + "_" + opts_name + "_result.txt"
   log_file = os.path.join(instance_dir, log_file)
   result_file = os.path.join(instance_dir, result_file)
   log_file = os.path.join(output_dir, log_file)
   result_file = os.path.join(output_dir, result_file)
   return log_file, result_file

def log_valid(instance_file, algorithm, output_dir, opts_name):
   log_file, result_file = get_log_file_names(instance_file, algorithm, output_dir, opts_name)
   # if not (os.path.isfile(log_file) and os.path.isfile(result_file)) or (os.path.isfile(log_file) and os.path.getsize(log_file) == 0) or (os.path.isfile(result_file) and os.path.getsize(result_file) == 0):
   if not (os.path.isfile(log_file)) or (os.path.isfile(log_file) and os.path.getsize(log_file) == 0):
       return False
   else:
       return True

def remove_incomplete_log(instance_file, algorithm, output_dir, opts_name):
   # remove log and result file if only one of them is present or either of them is empty
   log_file, result_file = get_log_file_names(instance_file, algorithm, output_dir, opts_name)

   if not log_valid(instance_file, algorithm, output_dir, opts_name):
      print "removing " + log_file + " and " + result_file
      os.remove(log_file)
      os.remove(result_file)

def evaluate_instance(instance_file, algorithm, executable_dir, instance_dir, output_dir, solver_options, opts_name):
   full_algorithm_path = os.path.join(executable_dir, algorithm)
   full_instance_path = os.path.join(instance_dir, instance_file)
   log_file, result_file = get_log_file_names(instance_file, algorithm, output_dir, opts_name)

   create_dir_cmd = "mkdir -p " + os.path.dirname(log_file) 
   call(create_dir_cmd, shell=True)

   if not (os.path.isfile(log_file) or os.path.isfile(result_file)):
      print("optimize " + instance_file + " with " + algorithm)
      run_alg_cmd = "export OMP_NUM_THREADS=1; " + full_algorithm_path + " " + full_instance_path + " " + solver_options + " 2>&1 > " + log_file
      print run_alg_cmd
      call(run_alg_cmd, shell=True)
#   else:
#      print("skip optimization of " + instance_file + " with algorithm " + algorithm)

def run_experiments(instance_list, algorithms, executable_dir, instance_dir, output_dir, solver_options, opts_names):
   for instance in instance_list:
      for algorithm, options, opts_name in zip(algorithms, solver_options, opts_names):
         evaluate_instance(instance, algorithm, executable_dir, instance_dir, output_dir, options, opts_name)

# RUN GUROBI

def evaluate_instance_gurobi(instance_file, instance_dir, output_dir):
   full_instance_path = os.path.join(instance_dir, instance_file)
   algorithm = 'gurobi_cl'
   params = "Presolve=0 LogToConsole=1 DisplayInterval=1 Threads=1 Method=1 NodeMethod=1 TimeLimit=36000"
   params_name = 'dual-simplex'
   log_file, result_file = get_log_file_names(instance_file, algorithm, output_dir, params_name)
   create_dir_cmd = "mkdir -p " + os.path.dirname(log_file) 
   call(create_dir_cmd, shell=True)

   if not (os.path.isfile(log_file) or os.path.isfile(result_file)):
      print("optimize " + instance_file + " with " + algorithm)
      run_alg_cmd = algorithm + " " + params + " " + full_instance_path + " 2>&1 > " + log_file
      print run_alg_cmd
      call(run_alg_cmd, shell=True)

def run_experiments_gurobi(instance_list, instance_dir, output_dir):
   for instance in instance_list:
      evaluate_instance_gurobi(instance, instance_dir, output_dir)

### EVALUATION ###

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
      if line.startswith("initial "):
         m = re.match("initial lower bound = (-?\d+\.?\d*), time = (\d+\.?\d*) s", line)
         if m:
            lower_bound = float(m.group(1))
            time_elapsed = float(m.group(2))
            lower_bounds.append(time_series_element(lower_bound, time_elapsed))
      if line.startswith("iteration "):
         m = re.match("iteration (\d+): lower bound = (-?\d+\.?\d*), time = (\d+\.?\d*) s", line)
         if m:
            iteration = float(m.group(1))
            lower_bound = float(m.group(2))
            time_elapsed = float(m.group(3))
            lower_bounds.append(time_series_element(lower_bound, time_elapsed))
      if line.startswith("Primal solution value:"):
         m = re.match("Primal solution value: (-?\d+\.?\d*), time = (\d+\.?\d*) s", line)
         if m:
            upper_bound = float(m.group(1))
            time_elapsed = float(m.group(2))
            upper_bounds.append(time_series_element(upper_bound, time_elapsed))
      if line.startswith("No primal solution found"):
         m = re.match("No primal solution found., time = (\d+\.?\d*) s", line)
         if m:
            upper_bound = 0.0
            time_elapsed = float(m.group(1))
            upper_bounds.append(time_series_element(upper_bound, time_elapsed))

   timeseries_continuous(lower_bounds)
   timeseries_continuous(upper_bounds)
   return lower_bounds, upper_bounds

def timeseries_continuous(timeseries):
   for i in range(1, len(timeseries)):
      if float(timeseries[i-1].time) > float(timeseries[i].time):
         print  str(i) + ": " + timeseries[i-1].time + "; " + timeseries[i].time
      assert float(timeseries[i-1].time) <= float(timeseries[i].time)

def average_time_series(time_series):
   no_time_series = len(time_series)
   time_series_concatenate = []
   for i in range(0,no_time_series):
      for elem in time_series[i]:
         time_series_concatenate.append(time_series_concatenation_element (i, elem.value, elem.time))
   sorted_ts = sorted(time_series_concatenate, key=lambda x: float(x.time))
   timeseries_continuous(sorted_ts)

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

   timeseries_continuous(averaged_ts)
   return averaged_ts

def average_plot(method_name, line_color, input_files):
   lower_bounds = []
   # upper_bounds = []
   for input_file in input_files:
      lb, ub = read_log(input_file)
      lower_bounds.append(lb)
      # upper_bounds.append(ub)

   averaged_lower_bounds = average_time_series(lower_bounds)
   # averaged_upper_bounds = average_time_series(upper_bounds)

   x_axis = []
   y_axis = []
   for x in averaged_lower_bounds:
      x_axis.append(float(x.time))
      y_axis.append(float(x.value))
   plt.plot(x_axis,y_axis, label=method_name, c=line_color)
   plt.legend()

   # x_axis = []
   # y_axis = []
   # for x in averaged_upper_bounds:
   #    x_axis.append(float(x.time))
   #    y_axis.append(float(x.value))
   # plt.plot(x_axis,y_axis, c=line_color, dashes=[6, 2])
   # plt.legend()

   plt.ylabel('average lower bound')
   plt.xlabel('time (s)')
   plt.tight_layout()


def average_results(input_files):
   lower_bounds = []
   upper_bounds = []
   times = []
   get_values = lambda l: map(lambda x: float(x.value), l)
   get_times = lambda l: map(lambda x: float(x.time), l)
   for input_file in input_files:
      lb, ub = read_log(input_file)
      lower_bounds.append(max(get_values(lb)))
      if get_values(ub):
         upper_bounds.append(min(get_values(ub)))
      else:
         upper_bounds.append(0.0);
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

def create_plots(log_files, methods, opts_names, output_dir, dataset_name):
   file_path = output_dir + dataset_name
   if os.path.isfile(file_path + '.txt'):
      os.remove(file_path + '.txt')
   with open(file_path + '.txt', 'w') as f:
      print >> f, dataset_name
      print >> f, "Avg UB | Avg LB | Avg time"
   # sort log files by method
   for method, opts_name, line_color in zip(methods, opts_names, line_colors):
      method_opts_name = method + "_" + opts_name
      print method_opts_name
      input_files = []
      for log_file in log_files:
         if log_file.endswith(method_opts_name + "_log.txt"): 
            input_files.append(log_file)

      average_plot(method_opts_name, line_color, input_files)
      #plt.show()
      #if path.endswith('/'): path = path[:-1]
      #figure_name = path.split('/')[-1]
      avg_lb, avg_ub, avg_time = average_results(input_files)
      print "Avg UB | Avg LB | Avg time"
      print avg_ub, avg_lb, avg_time
      with open(file_path + '.txt', 'a') as f:
         print >> f, method_opts_name
         print >> f, avg_ub, avg_lb, avg_time

   plt.savefig( file_path + '.pdf') 
   plt.close()

def create_table(paths, methods):
   table_latex = compute_primal_dual_table(methods, paths)
   table_latex = '\documentclass[preview]{standalone}\n\\usepackage{multirow}\n\\begin{document}\n' + table_latex + '\end{document}'
   open('results_table.txt', 'w').write(table_latex)
   print table_latex
   #call('pdflatex results_table.txt', shell=True)

def evaluate_experiments(instance_list, algorithms, opts_names, output_dir, dataset_name):
   log_files = []
   #print instance_list
   for instance in instance_list:
      print instance
      for algorithm, opts_name in zip(algorithms, opts_names):
         log_files.append(get_log_file_names(instance, algorithm, output_dir, opts_name)[0])
         if not log_valid(instance, algorithm, output_dir, opts_name):
            print log_files[-1] + " not valid, exit evaluation"
            exit(-1)
   print dataset_name
   create_plots(log_files, algorithms, opts_names, output_dir, dataset_name)
