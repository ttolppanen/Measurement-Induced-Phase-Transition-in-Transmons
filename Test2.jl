using ParametersModule, MIPTM

struct unmutableStruct
	a
	b
	function unmutableStruct(x)
		new(x, 2*x)
	end
end
mutable struct mutableStruct
	c
	d
end

function f()
	z(x) = 69*x
	a = [x->2*x, x->4*x]
	display(z(2))
end

f()
