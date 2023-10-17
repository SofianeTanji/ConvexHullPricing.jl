using JuMP

function Compute_A_B(u, nb_gen)
	nb_intervals_gen = zeros(Int64, nb_gen);
	if nb_gen>1
		for g=1:nb_gen
			flag = false;
			for elem in u[g,:]
				if elem>0
					if flag==false
						flag = true;
						nb_intervals_gen[g]+=1;
					end
				elseif flag==true
					flag = false;
				end
			end
		end
	else
		flag = false;
		for elem in u
			if elem>0
				if flag==false
					flag = true;
					nb_intervals_gen[1]+=1;
				end
			elseif flag==true
				flag = false;
			end
		end
	end
	A = (-1)*ones(Int64, nb_gen, maximum(nb_intervals_gen));
	B = (-1)*ones(Int64, nb_gen, maximum(nb_intervals_gen));
	a = 0;
	b = 0;
	if nb_gen>1
		for g=1:nb_gen
			flag = false;
			count = 0
			index = 1;
			for elem in u[g,:]
				count+=1;
				#println("elem : ",elem)
				if elem>0
					if flag==false
						flag = true;
						#A[g,index] = count;
						a = count;
					end
				elseif flag==true
					B[g,index] = count-1;
					A[g,index] = a;
					index+=1;
					flag = false;
				end
			end
		end
	else
		flag = false;
		count = 0
		index = 1;
		for elem in u
			count+=1;
			#println("[Compute_A_B] elem : ",elem)
			if elem>0
				if flag==false
					flag = true;
					#A[g,index] = count;
					a = count;
				end
			elseif flag==true
				B[1,index] = count-1;
				A[1,index] = a;
				index+=1;
				flag = false;
			end
		end
	end
	return A,B,nb_intervals_gen
end


function A_B_In_Intervals(a,b,g,A,B,nb_intervals_gen)
	for i=1:nb_intervals_gen[g]
		if a==A[g,i] && b==B[g,i]
			return true
		end
	end
	return false;
end


function Update_A_B(u, nb_gen, A, B, nb_intervals_gen)
	a = 0;
	b = 0;
	nb_intervals_added = zeros(Int64, nb_gen);
	nb_intervals_gen_old = copy(nb_intervals_gen);
	if nb_gen>1
		for g=1:nb_gen
			flag = false;
			count = 0;
			for elem in u[g,:]
				count+=1;
				if elem>0
					if flag==false
						flag = true;
						# nb_intervals_gen[g]+=1;
						a = count;
					end
				elseif flag==true
					b = count-1;
					if A_B_In_Intervals(a,b,g, A, B, nb_intervals_gen_old)==false
						nb_intervals_gen[g]+=1;
						nb_intervals_added[g]+=1;
					end
					flag = false;
				end
			end
		end
	else
		flag = false;
		count = 0;
		for elem in u
			count+=1;
			if elem>0
				if flag==false
					flag = true;
					# nb_intervals_gen[g]+=1;
					a = count;
				end
			elseif flag==true
				b = count-1;
				if A_B_In_Intervals(a,b,1, A, B, nb_intervals_gen_old)==false
					nb_intervals_gen[1]+=1;
					nb_intervals_added[1]+=1;
				end
				flag = false;
			end
		end
	end
	a = 0;
	b = 0;
	A_new = (-1)*ones(Int64, nb_gen, maximum(nb_intervals_gen));
	B_new = (-1)*ones(Int64, nb_gen, maximum(nb_intervals_gen));
	A_added = (-1)*ones(Int64, nb_gen, maximum(nb_intervals_added));
	B_added = (-1)*ones(Int64, nb_gen, maximum(nb_intervals_added));
	(i,j) = size(A);
	A_new[1:i,1:j] = A;
	(i,j) = size(B);
	B_new[1:i,1:j] = B;
	if nb_gen>1
		for g=1:nb_gen
			flag = false;
			index = nb_intervals_gen_old[g]+1;
			index_added = 1;
			count = 0;
			for elem in u[g,:]
				count+=1;
				if elem>0
					if flag==false
						flag = true;
						a = count;
					end
				elseif flag==true
					b = count-1;
					if A_B_In_Intervals(a,b,g, A, B, nb_intervals_gen_old)==false
						A_new[g,index] = a;
						B_new[g,index] = b;
						A_added[g,index_added] = a;
						B_added[g,index_added] = b;
						index+=1;index_added+=1;
					end
					flag = false;
				end
			end
		end
	else
		#println("One generator");
		flag = false;
		index = nb_intervals_gen_old[1]+1;
		index_added = 1;
		count = 0;
		for elem in u
			count+=1;
			if elem>0
				if flag==false
					flag = true;
					a = count;
				end
			elseif flag==true
				b = count-1;
				if A_B_In_Intervals(a,b,1, A, B, nb_intervals_gen_old)==false
					A_new[1,index] = a;
					B_new[1,index] = b;
					A_added[1,index_added] = a;
					B_added[1,index_added] = b;
					index+=1;index_added+=1;
				end
				flag = false;
			end
		end
	end
	return A_new,B_new,nb_intervals_gen, A_added, B_added, nb_intervals_added
end

function ColumnAndRowGeneration(instance, niter, eps)
end