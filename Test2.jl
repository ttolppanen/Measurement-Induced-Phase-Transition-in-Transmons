using RecursiveArrayTools

function test(i)
	2*i
end

function arrayForEveryThread()
	a = []
	for _ in 1:Threads.nthreads()
		push!(a, [])
	end
	a
end

function f(max)
	a = arrayForEveryThread()
	Threads.@threads for i = 1:max
		push!(a[Threads.threadid()], test(i))
	end
	a = reduce(vcat, a)
	sum(a)
end

function fb(max)
	a = []
	Threads.@threads for i = 1:max
		push!(a, test(i))
	end
	sum(a)
end

display(f(100))
