import math
import numpy

# parser = argparse.ArgumentParser(description="Naive Bayes model using a training set and predict the tissues of a testing set with its accuracy.")

# parser.add_argument('-i', '--input', 
# 					dest = "infile", 
# 					action = "store",
# 					default = None,
# 					help = "Input file")

# parser.add_argument('-o', '--output', 
# 					dest = "outfile", 
# 					action = "store",
# 					default = None,
# 					help = "Output file")

# parser.add_argument('-v', '--verbose', 
# 					dest = "verbose", 
# 					action = "store_true",
# 					default = False,
# 					help = "Print log in stderr")

# options = parser.parse_args()

def NaiveBayes(training_file, testing_file, psi_file):
	# Variables
	training_dic = {} #key (tissue) -> value (sampleid)
	testing_dic = {} #key (tissue) -> value (sampleid)

	training_event_dic = {} #key (event) -> value (values_up, values_down lists)
	training = {} #key (tissue=key of training_dic) -> value (training_event_dic=event:values_up, values_down)

	samples_testing_list = []


	################### Reading training_file ###################
	fd = open(training_file, 'r')
	for line in fd:
		fields = line.rstrip('\n').split('\t')
		tissue = fields[2]
		sampleid = fields[0]
		training_dic.setdefault(tissue, []).append(sampleid)
	# return(training_dic)
	fd.close()

	################### Reading testing_file ###################
	fd = open(training_file, 'r')
	for line in fd:
		fields = line.rstrip('\n').split('\t')
		tissue = fields[2]
		sampleid = fields[0]
		samples_testing_list.append(sampleid)
		testing_dic.setdefault(tissue, []).append(sampleid)
	fd.close()

	################### Reading psi_file ###################
	fd = open(psi_file, 'r')
	line_sample = fd.readline() #read first line (which has the samples id)
	samples = line_sample.split('\t')
	# return len(samples) #1477

	for key,value in training_dic.items():
		training_event_dic = {}
		fd.seek(0)
		next(fd)
		# return(training_event_dic)

		for line in fd:
			# return(line)
			fields = line.rstrip('\n').split('\t')
			event = fields[0]
			values = fields[1:]
			values_up = []
			values_down = []

			i=0
			for sample in value: #value of training_dic (sampleid)
				i = samples.index(sample)
				if values[i] !="NA":
					if float(values[i]) == 1: 
						values_up.append(samples[i])
					elif float(values[i]) == 0: 
						values_down.append(samples[i])
			training_event_dic[event] = [["up", values_up], ["down", values_down]]
		training[key] = training_event_dic #key of training_dic (tissue)
		# return(len(training_event_dic.keys())) #20373 events
		return(training)

	for key,value in testing_dic.items():
		testing_event_dic = {}

result = NaiveBayes("training_set.description", "testing_set.description", "psi.formatted.test")
print(result)



