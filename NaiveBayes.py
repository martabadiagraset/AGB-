import math
import numpy
import os

def NaiveBayes(training_file, testing_file, psi_file):

#####################################################################
######### Preparing Training and Testing dictinaries ################
#####################################################################

	########################### Variables ###########################
	tissues_set = set()
	tissues = []
	training_dic = {} #key (tissue) -> value (sampleid)
	testing_dic = {} #key (tissue) -> value (sampleid)

	training_event_dic = {} #key (event) -> value (values_up, values_down lists)
	training = {} #key (tissue=key of training_dic) -> value (training_event_dic=event:values_up, values_down)

	values = []
	list_samples = [] #all samples id
	samples_testing_list = []
	samples_testing = []

	testing_sample_dic = {}
	testing_tissue = {}
	testing = {}

	######################## Reading training_file ####################
	with open(training_file, 'r') as fd:
		for line in fd:
			fields = line.rstrip('\n').split('\t')
			tissue = fields[2]
			sampleid = fields[0]
			tissues_set.add(tissue)
			tissues = list(tissues_set)
			training_dic.setdefault(tissue, []).append(sampleid)
		# return(training_dic)

	######################### Reading testing_file ######################
	with open(testing_file, 'r') as fd:
		for line in fd:
			fields = line.rstrip('\n').split('\t')
			tissue = fields[2]
			sampleid = fields[0]
			samples_testing_list.append(sampleid)
			testing_dic.setdefault(tissue, []).append(sampleid)

	########################## Reading psi_file #########################
	fd = open(psi_file, 'r')
	line_sample = fd.readline() #read first line (which has the samples id)
	samples = line_sample.split('\t')
	# return len(samples) #1477

	for sample in samples:
		list_samples.append(sample.rstrip('\n')) 

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
				i = list_samples.index(sample)
				if values[i] !="NA":
					if float(values[i]) == 1: 
						values_up.append(list_samples[i])
					elif float(values[i]) == 0: 
						values_down.append(list_samples[i])
			training_event_dic[event] = [["up", values_up], ["down", values_down]]
		training[key] = training_event_dic #key of training_dic (tissue)
		# return(len(training_event_dic.keys())) #20373 events
		# return(training)

	for key,value in testing_dic.items():
		testing_sample_dic = {}
		samples_testing = []

		for sample in value: #value of testing_dic (sampleid)
			i=0
			samples_testing.append(sample)
			i = list_samples.index(sample)
			SE_up = []
			SE_down = []

			with open(psi_file, 'r') as fd:
				next(fd)
				for line in fd:
					fields = line.rstrip('\n').split('\t')
					event = fields[0]
					values = fields[1:]
					# return(event,values)
					if values[i] !="NA":
						if float(values[i]) == 1: 
							SE_up.append(event)
						elif float(values[i]) == 0: 
							SE_down.append(event)
				testing_sample_dic[sample] = [["up", SE_up], ["down", SE_down]] #key of testing_sample_dic (sample)
		testing_tissue[key] = samples_testing

		for key, values in testing_sample_dic.items():
			testing[key] = values

#####################################################################
########## Likelihood: prior probabilities -> P(SE|tissue) ##########
#####################################################################

	samples = 20
	value_total = {} #tissue -> SE -> up: count, down: count
	prob_training = {} #tissue -> SE -> up: prob, down: prob

	for tissue, value in training.items(): #for tissue
		value_total_tissue = {}
		prob_training_tissue = {}

		for SE, value2 in value.items(): #for SE
			value_SE = []
			prob_SE = []

			for updown in value2: #for up, down
				value_SE.append([updown[0], len(updown[1])]) #total number of up and down events (attribute=up and down)
				prob_SE.append([updown[0], (len(updown[1])+1)/(samples+2)])

			value_total_tissue[SE] = value_SE #Nº of samples for each SE(up/down) in each tissue (key=event / value=up and down counts -> nº of samples)
			prob_training_tissue[SE] = prob_SE #Likelihood (key=event / value=up and down probabilities)

		value_total[tissue] = value_total_tissue
		prob_training[tissue] = prob_training_tissue
		# return tissue, SE, value_SE, prob_SE
	# return value_total
	# return prob_training

#####################################################################
####### Posterior probabilities: probability of each tissue #########
####### given the attribute (SE) -> P(tissue|SE) P(SE) ##############
#####################################################################

	PSE_dic = {}
	Ptissue_SE = {}

	# tissues=["Liver", "Heart", "Nerve", "Lung", "Muscle"]

	tissues_total = [value_total[tissues[0]], value_total[tissues[1]], value_total[tissues[2]], value_total[tissues[3]], value_total[tissues[4]]] #push all content (value= event:up[] down[] values) of value_total dictionary for all tissues (as tissues is the key of value_total)
	#tissues_total = tissue content (SE, up and down values)

	# return value_total[tissues[0]]
	# return tissues_total

	tissues_training = [prob_training[tissues[0]], prob_training[tissues[1]], prob_training[tissues[2]], prob_training[tissues[3]], prob_training[tissues[4]]] #push all content (value= event:up[] down[] probabilities) of prob_training dictionary for all tissues (as tissues is the key of prob_training)
	# return tissues_training
	#tissues_training = tissue content (SE, up and down probabilities)

	# SEs in all tissues (unique)
	SE_uniques = set(value_total[tissues[0]].keys()).intersection(set(value_total[tissues[1]].keys()), set(value_total[tissues[2]].keys()), set(value_total[tissues[3]].keys()), set(value_total[tissues[4]].keys()))
	# return value_total[tissues[0]].keys()
	# return SE_unique

	for SE in SE_uniques:
		up_prob = 0.0
		down_prob = 0.0
		list_prob = []

		for tissue in tissues_total:
			if SE in tissue:
				up_prob+=tissue[SE][0][1]
				down_prob+=tissue[SE][1][1]
		#P(SE)
		PSE_dic[SE]=[["up",((up_prob+1)/(100+2))], ["down", ((down_prob+1)/(100+2))]] 

		for tissue in tissues_total:
			list_prob.append([["up", ((tissue[SE][0][1]+1)/(up_prob+5))], ["down", ((tissue[SE][1][1]+1)/(down_prob+5))]]) #list_prob = list of SE up and down probabilities in all tissues
		
		#P(tissue|SE): Posterior probabilities -> key=SE / value=list_prob (list of all up and down probabilities of each SE in all tissues) 	
		Ptissue_SE.setdefault(SE,list_prob)
		# return Ptissue_SE 

#####################################################################
########## Mutual Information (Information Gain, MI or IG) ##########
#####################################################################
	
	n_tissues = len(tissues)
	
	# Entropies (H). H(tissue) and H(tissue|SE)
	Htissue = n_tissues*-1*((1/n_tissues)*math.log2(1/n_tissues)) #H(tissue)
	Htissue_SE = {} #H(tissue|SE) = key(SE), value (H)

	# Mutual information (MI). MI(tissue,SE)=H(tissue)-H(tissue|SE)
	MItissue_SE = {} #MI(tissue, SE) = key(SE), value (MI)

	for SE, probs_tissue in Ptissue_SE.items():
		H_up = 0
		H_down = 0

		for tissue in probs_tissue:
			H_up+=tissue[0][1]*math.log2(tissue[0][1])
			H_down+=tissue[1][1]*math.log2(tissue[1][1])

		Htissue_SE[SE] = ((-1*PSE_dic[SE][0][1]*H_up) + (-1*PSE_dic[SE][1][1]*H_down))
		MItissue_SE[SE] =  Htissue - Htissue_SE[SE]

		#return MItissue_SE
	# return MItissue_SE
	# return PSE_dic
	
	### Outputs with ALL SE and its corresponent H and MI values -> merge them in R by SE column ###
	# Output file with: all SE and its corresponent H value (SE sorted by H values)
	entropy_output = ("output/Entropy.txt")
	if os.path.exists(entropy_output):
		print ("Entropy.txt has been already created!")
	else:
		for SE in (sorted(Htissue_SE, key=Htissue_SE.get, reverse=True)):
			with open (entropy_output, "a") as out:
				out.write("%s\t%d\n"%(SE, Htissue_SE[SE]))

	# Output file with: all MI values sorted
	MI_output = ("output/Mututal_information.txt")
	if os.path.exists(MI_output):
		print ("Mututal_information.txt has been already created!")
	else:
		for MI in (sorted(MItissue_SE.values(), reverse=True)):
			with open(MI_output, "a") as out:
				out.write("%s\n"%(MI))

	# Output file with: all SE and its corresponent MI value (SE sorted by MI values)
	MI_SE_output = ("output/Mututal_information_SE.txt")
	if os.path.exists(MI_SE_output):
		print ("Mututal_information_SE.txt has been already created!")
	else:
		for SE in (sorted(MItissue_SE, key=MItissue_SE.get, reverse=True)):
			with open(MI_SE_output, "a") as out:
				out.write("%s\t%s\n"%(SE, MItissue_SE[SE]))

	### Selecting SE depending on their MI value (define a cut-off), from MItissue_SE ###
	MI_values = numpy.array(list(map(float, MItissue_SE.values())))
	MI_cutoff = numpy.percentile(MI_values, 75) #to use percentile you need a numpy.array!
	# return MI_cutoff

	SE_selected = []
	for SE,MI in MItissue_SE.items():
		if MI >= MI_cutoff:
			SE_selected.append(SE)
	print ("Total number of SE selected: %d (MI cut-off %.4f)"%(len(SE_selected), MI_cutoff))

#####################################################################
########## Naive Bayes (using the testing samples) ##################
#####################################################################

	NB_results = []
	TP = 0
	FP = 0
	# return tissues_training
	# return testing
	for sample_test in samples_testing_list:
		sample_test_prob = []

		for tissue in tissues_training:
			tissue_prob = 0.0

			for SE in SE_selected:
				# return tissues_training
				# return tissue[SE]
				if SE in testing[sample_test][0][1]: #look if in the sample_test, the SE_selected is in list of up events
					tissue_prob+=math.log2(tissue[SE][0][1])

				elif SE in testing[sample_test][1][1]: #look if in the sample_test, the SE_selected is in list of down events
					tissue_prob+=math.log2(tissue[SE][1][1])
				# return tissue_prob
			print(tissue[SE])
			print(tissue_prob)
			sample_test_prob.append(tissue_prob*1/n_tissues) #list of the probabilities of the sample_test to correspond to the n tissues
		print (sample_test_prob)

		print("%s done"%sample_test)

		score = max(sample_test_prob) #highest probability of the sample (corresponding to one of the n tissues)
		
		# Find the predicted_tissue (prediction)
		score_index = sample_test_prob.index(score)
		predicted_tissue = tissues[score_index]

		# Find the real_tissue (label)
		for tissue in testing_tissue:
			if sample_test in testing_tissue[tissue]:
				label_tissue = tissue

		NB_results.append((score, predicted_tissue, label_tissue, sample_test))

	# return NB_results
	NB_results.sort(reverse=True) #sort by score

	# Output file with: NB prediction results
	NB_results_output = ("output/NB_results.txt")
	if os.path.exists(NB_results_output):
		print ("NB_results.txt has been already created!")
	else:
		with open(NB_results_output, "a") as out_results:
			out_results.write("Total number of SE selected: %d.\nMutual Information cut-off %.4f\n\nTesting sample\tScore\tPrediction\tLabel\n"%(len(SE_selected), MI_cutoff))
			for i in NB_results:
				results_table = str("%s\t%.3f\t%s\t%s\n"%(i[3], i[0], i[1], i[2]))
				out_results.write(results_table)

#####################################################################
#################### Accuracy of the models #########################
#####################################################################
	
	accuracy = {}
	score_cutoff_list = [-400, -600, -800, -1000, -1200, -1400, -2000]

	for tissue in tissues:
		for score in score_cutoff_list:
			accuracy[(score, tissue)] = [["TP", 0], ["FP", 0], ["FN", 0], ["TN", 0], ["TPR", 0], ["FPR", 0], ["FDR", 0]]

	# return accuracy[(score, tissue)]

	for item in NB_results: #items=line
		for score in score_cutoff_list:
			if item[0] > score:
				for tissue_label in tissues:
					if item[1]==tissue_label and item[2]==tissue_label: #TP (OK)
						accuracy[(score, tissue_label)][0][1]+=1

					elif item[1]==tissue_label and item[2]!=tissue_label: #FP 
						accuracy[(score, tissue_label)][1][1]+=1

					elif item[2]==tissue_label and item[1]!=tissue_label: #FN
						accuracy[(score, tissue_label)][2][1]+=1

					else:
						accuracy[(score, tissue_label)][3][1]+=1 #TN (OK)

	for score_tissue in accuracy: #score_tissue is key of accuracy = (score, tissue)
		# True Positive Rate (TPR) -TP and FN-
		if accuracy[score_tissue][0][1]!=0 or accuracy[score_tissue][2][1]!=0:
			TP = accuracy[score_tissue][0][1]
			FN = accuracy[score_tissue][2][1]
			TPR = TP/(TP+FN)
			accuracy[score_tissue][4][1] = TPR

		# False Positive Rate (FPR) -TN and FP-
		if accuracy[score_tissue][3][1]!=0 or accuracy[score_tissue][1][1]!=0:
			TN = accuracy[score_tissue][3][1]
			FP = accuracy[score_tissue][1][1]
			FPR = FP/(FP+TN)
			accuracy[score_tissue][5][1] = FPR

		# False Discovery Rate (FDR) -FP and TP-
		if accuracy[score_tissue][1][1]!=0 or accuracy[score_tissue][0][1]!=0:
			FP = accuracy[score_tissue][1][1]
			TP = accuracy[score_tissue][0][1]
			FDR = FP/(FP+TP)
			accuracy[score_tissue][6][1] = FDR

	Accuracy_output = ("output/Accuracy.out")
	if os.path.exists(Accuracy_output):
		print ("Accuracy.out has been already created!")
	else:
		with open(Accuracy_output, "a") as out_accuracy:
			out_accuracy.write("Score\tTissue\tTP\tFP\tFN\tTN\tTPR\tFPR\tFDR\n")
			for key, value in accuracy.items():
				out_accuracy.write("%s\t%s\t%s\t%s\t%s\t%s\t%.3f\t%.3f\t%.3f\n"%(key[0],key[1], value[0][1], value[1][1],value[2][1],value[3][1],value[4][1],value[5][1], value[6][1]))

# parse_input("training_set", "testing_set", "psi.formatted.test")

result = NaiveBayes("test_training_set", "test_testing", "psi.formatted")





