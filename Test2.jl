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
	sp = BoseHubbardParameters(L=3, N=3, U=0.14)
	display(sp.isThereDisorderInW)
end

f()
