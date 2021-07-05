
function f()
	a = []
	b = []
	for i in 1:Threads.nthreads()
		push!(a, [])
		push!(b, [])
	end
	Threads.@threads for i in 1:100
		c = i*3
		push!(a[Threads.threadid()], c)
		push!(b[Threads.threadid()], c^2)
	end
	a = reduce(vcat, a)
	b = reduce(vcat, b)
	print(sum(b))
end
f()
