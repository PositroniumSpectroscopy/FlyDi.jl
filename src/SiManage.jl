#-----Functions to manage simulations-----
function MakDir()
	sid = Dates.format(Dates.now(), "yyyymmdd_HHMMSS")
	newdir = joinpath("Z:\\Decel - Simulation", sid)
	mkdir(newdir)
	newdir = joinpath("C:\\Users\\ucapesh\\Desktop\\OneDrive - University College London\\Simulation Data",sid)
	mkdir(newdir)
	return(sid)
end

function DFSave(sid::String, fname::String, data::Array{DataFrame,1})
	path = joinpath("Z:\\Decel - Simulation", sid)
	DatOut = joinpath(path, fname)
	Feather.write(DatOut, data)
	
	path = joinpath("C:\\Users\\ucapesh\\Desktop\\OneDrive - University College London\\Simulation Data",sid)
	DatOut = joinpath(path, fname)
	Feather.write(DatOut, data)
end

function DFOpen(sid::String, fname::String)
	path = joinpath("Z:\\Decel - Simulation", sid)
	DatIn = joinpath(path, fname)
	data = Feather.read(DatIn)::DataFrame
	return data
end
