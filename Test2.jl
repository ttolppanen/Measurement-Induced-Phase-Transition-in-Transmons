using Plots

function f()
	H = [1 2; 3 4]
	res = []
	for _ in 1:Threads.nthreads()
		push!(res, 0)
	end
	Threads.@threads for _ in 1:1000
		res[Threads.threadid()] += [1, 4]' * H * [1, 2]
	end
	display(typeof(H))
	display(sum(res))
end

f()
