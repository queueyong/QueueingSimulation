function time_axis(file_array::Matrix{Any})
 temp_array = file_array[1,:]
 x = Float64[]
 i = 1
 while typeof(temp_array[i]) != SubString{String}
  push!(x, temp_array[i])
  i += 1
 end
 return x
end

function time_axis(file_array::Matrix{Float64})
 temp_array = file_array[1,:]
 return temp_array
end

function save_mean_with_ci(file_array::Any, stream::IOStream)
 x = time_axis(file_array)
 avg = Float64[]
 upper_ci = Float64[]
 lower_ci = Float64[]
 replication = size(file_array)[1]-1
 for i in 1:length(x)
  values = convert(Array{Float64,1}, file_array[2:end,i])
  σ = std(values)
  push!(avg, sum(values)/replication)
  push!(upper_ci, avg[i]+1.96*(σ/sqrt(replication)))
  push!(lower_ci, avg[i]-1.96*(σ/sqrt(replication)))
 end
 writedlm(stream, transpose(x)) # time axis
 writedlm(stream, transpose(avg)) # average
 writedlm(stream, transpose(upper_ci)) # upper ci
 writedlm(stream, transpose(lower_ci)) # lower ci
end

function save_tail_prob(file_array::Any, sum_tail::IOStream, δ_set::Vector{Float64})
 x = time_axis(file_array)
 writedlm(sum_tail, transpose(x))
 for δ in δ_set
  prob = Float64[]
  for j in 1:length(x)
   sum = 0
   for i in 2:(size(file_array)[1]-1)
    if file_array[i,j] > δ
      sum += 1
    end
   end
   push!(prob, sum/(size(file_array)[1]-1) )
  end
  writedlm(sum_tail, transpose(prob))
 end
end

queue = "TVGG1PS"
#control_set = ("RM", "SR", "PSA")
#control_set = ("RM", "SR", "PSA", "DM")
control_set = ("RM",)
#control_set = ("DM",)
#param_set = (0.8, 0.2, 1.0)
#param_set = (0.8, 0.2, 1.0, 0.25)
param_set = (0.8,)
#param_set = (0.25,)
#param_letter_set = ("ρ", "ξ", "ζ")
#param_letter_set = ("ρ", "ξ", "ζ", "δ")
param_letter_set = ("ρ",)
#param_letter_set = ("δ",)

distribution_set = [("Exponential","Exponential"), ("Hyperexponential","Hyperexponential"), ("Hyperexponential","Erlang"), ("Erlang","Erlang")]
gamma_set = (0.001, 0.01, 0.1, 1.0 , 10.0)
T_set = (20000.0, 2000.0, 2000.0, 2000.0, 2000.0)
RAW_LOG_FILE_DIR = "C:/Users/YOC/Google 드라이브/YOC-Data/Simulation Log/$queue"
# RAW_LOG_FILE_DIR = "D:/Julia/$queue/Log files"
SAVE_DIR = "C:/Users/YOC/Google 드라이브/YOC-Data/Summarized Simulation Log/$queue"
rep = 10000
δ_set = [1.0, 3.0, 5.0, 7.0, 9.0, 11.0, 13.0]

# mean num & mean wait
for c in 1:length(control_set)
 for d in 1:length(distribution_set)
  for γ in 1:length(gamma_set)
   temp_string = string(RAW_LOG_FILE_DIR,"/",control_set[c],"/",string(distribution_set[d][1],", ",distribution_set[d][2]),"/")
   file_name_1 = "num in queue ($queue, gamma $(gamma_set[γ]), control $(control_set[c]), param $(param_set[c]), arrival $(distribution_set[d][1]), service $(distribution_set[d][2]), time $(T_set[γ]), rep $rep).txt"
   file_name_2 = "virtual sojourn time ($queue, gamma $(gamma_set[γ]), control $(control_set[c]), param $(param_set[c]), arrival $(distribution_set[d][1]), service $(distribution_set[d][2]), time $(T_set[γ]), rep $rep).txt"
   file_path_1 = string(temp_string,file_name_1)
   file_path_2 = string(temp_string,file_name_2)
   num_in_queue_file = open(file_path_1)
   virtual_sojourn_time_file = open(file_path_2)
   num_in_queue_file_array = readdlm(num_in_queue_file)
   virtual_sojourn_time_file_array = readdlm(virtual_sojourn_time_file)
   close(num_in_queue_file)
   close(virtual_sojourn_time_file)

   temp_string = string(SAVE_DIR,"/",control_set[c],"/",string(distribution_set[d][1],", ",distribution_set[d][2]),"/")
   file_name_3 = "(sum) num in queue ($queue, gamma $(gamma_set[γ]), control $(control_set[c]), arrival $(distribution_set[d][1]), service $(distribution_set[d][2]), time $(T_set[γ])).txt"
   file_name_4 = "(sum) virtual sojourn time ($queue, gamma $(gamma_set[γ]), control $(control_set[c]), arrival $(distribution_set[d][1]), service $(distribution_set[d][2]), time $(T_set[γ])).txt"
   sum_num = open(string(temp_string,file_name_3) ,"w")
   sum_sojourn = open(string(temp_string,file_name_4) , "w")
   save_mean_with_ci(num_in_queue_file_array, sum_num)
   save_mean_with_ci(virtual_sojourn_time_file_array, sum_sojourn)
   close(sum_num)
   close(sum_sojourn)

   file_name_5 = "(sum) tail probability ($queue, gamma $(gamma_set[γ]), control $(control_set[c]), arrival $(distribution_set[d][1]), service $(distribution_set[d][2]), time $(T_set[γ])).txt"
   sum_tail = open(string(temp_string,file_name_5) , "w")
   writedlm(sum_tail, transpose(δ_set))
   save_tail_prob(virtual_sojourn_time_file_array, sum_tail, δ_set)
   close(sum_tail)
  end
 end
end

# tail prob
for c in 1:length(control_set)
 for d in 1:length(distribution_set)
  for γ in 1:length(gamma_set)
   temp_string = string(RAW_LOG_FILE_DIR,"/",control_set[c],"/",string(distribution_set[d][1],", ",distribution_set[d][2]),"/")
   file_name_2 = "virtual sojourn time ($queue, gamma $(gamma_set[γ]), control $(control_set[c]), param $(param_set[c]), arrival $(distribution_set[d][1]), service $(distribution_set[d][2]), time $(T_set[γ]), rep $rep).txt"
   file_path_2 = string(temp_string,file_name_2)
   virtual_sojourn_time_file = open(file_path_2)
   virtual_sojourn_time_file_array = readdlm(virtual_sojourn_time_file)
   close(virtual_sojourn_time_file)
   temp_string = string(SAVE_DIR,"/",control_set[c],"/",string(distribution_set[d][1],", ",distribution_set[d][2]),"/")
   file_name_5 = "(sum) tail probability ($queue, gamma $(gamma_set[γ]), control $(control_set[c]), arrival $(distribution_set[d][1]), service $(distribution_set[d][2]), time $(T_set[γ])).txt"
   sum_tail = open(string(temp_string,file_name_5) , "w")
   writedlm(sum_tail, transpose(δ_set))
   save_tail_prob(virtual_sojourn_time_file_array, sum_tail, δ_set)
   close(sum_tail)
  end
 end
end
