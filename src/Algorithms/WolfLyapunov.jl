using DynamicalSystems, LinearAlgebra

export WolfAlgorithm
export solve

#struct to hold data needed for wolf algorithm and to allow for dispatch when solving 
struct WolfAlgorithm <: LCEAlgorithm
    ires
    maxbox
    dt 
    evolve 
    dismin 
    dismax 
    thmax 
end

function WolfAlgorithm(;ires = 10, maxbox = 6000, dt = 0.01, evolve = 20, dismin = 0.001, dismax = 0.3, thmax = 30)
    WolfAlgorithm(ires,maxbox,dt,evolve,dismin,dismax,thmax)
end
#struct to hold basegen data base
struct Basegen_db
    ndim 
    ires 
    tau 
    datcnt 
    boxcnt 
    datmax 
    datmin 
    boxlen 
    datptr 
    nxtbox 
    wherear
    nxtdat
    data
end

#this is a direct translation of Matlab code for the Wolf algorithm to Julia code, 
#from https://www.mathworks.com/matlabcentral/fileexchange/48084-wolf-lyapunov-exponent-estimation-from-a-time-series


# search function
function  search(iflag, ndim, ires, datmin,
    boxlen, nxtbox, wherear, datptr, nxtdat, data, delay, oldpnt, newpnt,
    datuse, dismin, dismax, thmax, evolve)
# Searches for the most viable point for fet.m
# Taehyeun Park; The Cooper Union; EE'15

target = zeros(1,ndim)
oldcrd = zeros(1,ndim)
zewcrd = zeros(1,ndim)

oldcrd[1:ndim] = data[Int.(oldpnt .+ delay)]
zewcrd[1:ndim] = data[Int.(newpnt .+ delay)]
igcrds = floor.((oldcrd .- datmin)./boxlen)
oldist = sqrt(sum((oldcrd .- zewcrd).^2))

irange = round(dismin/boxlen)
if irange .== 0
    irange = 1
end

thbest = thmax
bstdis = dismax
bstpnt = 0

goto30 = 1
while goto30 .== 1
    goto30 = 0
    for icnt = 0:((2*irange+1)^ndim)-1
        goto140 = 0
        icounter = icnt
        for ii = 1:ndim
            ipower = (2*irange+1)^(ndim-ii)
            ioff = floor(icounter./ipower)
            icounter = icounter - ioff*ipower
            target[ii] = igcrds[ii] - irange + ioff

            if target[ii] .< 0
                goto140 = 1
                break
            end
            if target[ii] .> ires-1
                goto140 = 1
                break
            end
        end
        
        if goto140 .== 1
            continue
        end
        
        if irange != 1
            iskip = 1
            for ii = 1:ndim
                if abs(round(target[ii] - igcrds[ii])) .== irange
                    iskip = 0
                end
            end
            if iskip .== 1
                continue
            end
        end
        
        runner = 1
        for ii = 1:ndim
            goto80 = 0
            goto70 = 1
            while goto70 .== 1
                goto70 = 0
                if wherear[Int(runner),ii] .== target[ii]
                    goto80 = 1
                    break
                end
                runner = nxtbox[Int(runner), ii]
                if runner != 0
                    goto70 = 1
                end
            end
            
            if goto80 .== 1
                continue
            end
            goto140 = 1
            break
        end
        
        if goto140 .== 1
            continue
        end
        
        if runner .== 0
            continue
        end
        runner = datptr[Int(runner)]
        if runner .== 0
            continue
        end
        goto90 = 1
        while goto90 .== 1
            goto90 = 0
            while true
                if abs(round(runner - oldpnt)) .< evolve
                    break
                end
                if abs(round(runner - datuse)) .< (2*evolve)
                    break
                end
                
                bstcrd = data[Int.(runner .+ delay)]
                
                abc1 = oldcrd[1:ndim] - bstcrd[1:ndim]
                abc2 = oldcrd[1:ndim] - zewcrd[1:ndim]
                tdist = sum(abc1.*abc1)
                tdist = sqrt(tdist)
                dot = sum(abc1.*abc2)

                if tdist .< dismin
                    break
                end
                if tdist >= bstdis
                    break
                end
                if tdist .== 0
                    break
                end
                goto120 = 0
                if iflag .== 0
                    goto120 = 1
                end
                if goto120 .== 0
                    ctheta = min(abs(dot/(tdist*oldist)),1)
                    theta = 57.3*acos(ctheta)
                    if theta >= thbest
                        break
                    end
                    thbest = theta
                end
                bstdis = tdist
                bstpnt = runner
                break
            end
            runner = nxtdat[Int(runner)]

            if runner != 0
                goto90 = 1
            end
        end
    end
    irange = irange + 1
    if irange <= (0.5 + round((dismax/boxlen)))
        goto30 = 1
        continue
    end
    
end
return bstpnt,bstdis,thbest
end




#basegen
function basgen(fname::String, tau, ndim, ires, datcnt, maxbox)
    # Database generator for fet.m function
    # Taehyeun Park; The Cooper Union; EE'15
    
    x = readlines(fname)
    #data = zeros(1,datcnt)
    data = parse.(Float64,x)
    trck = 1
    start = 1
    fin = 0
    
    #for ii = 1:length(x)
    #    #if cmp(x[ii], char(32)) || cmp(x[ii], Char(13)) || cmp(x[ii], Char(10)) || cmp(x[ii], Char(26))
    #    if x[ii] == Char(32) || x[ii] == Char(13) || x[ii] == Char(10) || x[ii] == Char(26)
    #        if fin >= start()
    #            data[trck] = parse.(Float64,x[start:fin])
    #            trck = trck + 1
    #            if trck .> 8*floor(datcnt/8)
    #                break
    #            end
    #        end
    #        start = ii + 1
    #    else
    #        fin = ii
    #    end
    #end




    
    delay = collect(0:tau:(ndim-1)*tau)
    
    nxtbox = zeros(maxbox, ndim)
    wherear = zeros(maxbox, ndim)
    datptr = zeros(1,maxbox)
    nxtdat = zeros(1,datcnt)
    
    datmin = minimum(data)
    datmax = maximum(data)
    
    datmin = datmin - 0.01*(datmax - datmin)
    datmax = datmax + 0.01*(datmax - datmin)
    boxlen = (datmax - datmin)/ires
    
    boxcnt = 1
    
    for ii = 1:(datcnt-(ndim-1)*tau)
        target = floor.((data[ii .+ delay] .- datmin)/boxlen)
        runner = 1
        chaser = 0
        
        jj = 1
        while jj <= ndim
            tmp = wherear[Int(runner),jj]-target[jj]
            if tmp < 0
                chaser = runner
                runner = Int(nxtbox[runner,jj])
                if runner != 0
                    continue
                end
            end
            if tmp != 0
               boxcnt = boxcnt + 1
               
               if boxcnt .== maxbox
                   error("Grid overflow, increase number of box count")
               end
               
               for kk = 1:ndim
                   wherear[boxcnt,kk] = wherear[chaser,kk]
               end
               wherear[boxcnt,jj] = target[jj]
               nxtbox[chaser,jj] = boxcnt
               nxtbox[boxcnt,jj] = runner
               runner = Int(boxcnt)
            end
            jj = jj + 1
        end
        nxtdat[ii] = datptr[runner]
        datptr[runner] = ii
    end
    
    used = 0
    for ii = 1:boxcnt
        if datptr[ii] != 0
            used = used + 1
        end
    end
    db = Basegen_db(ndim,ires,tau,datcnt,boxcnt,datmax,datmin,boxlen,datptr[1:boxcnt],nxtbox[1:boxcnt,1:ndim],wherear[1:boxcnt,1:ndim],nxtdat[1:datcnt],data)
    #db.ndim = ndim
    #db.ires = ires
    #db.tau = tau
    #db.datcnt = datcnt
    #db.boxcnt = boxcnt
    #db.datmax = datmax
    #db.datmin = datmin
    #db.boxlen = boxlen
    
    #db.datptr = datptr[1:boxcnt]
    #db.nxtbox = nxtbox[1:boxcnt, 1:ndim]
    #db.wherear = wherear[1:boxcnt, 1:ndim]
    #db.nxtdat = nxtdat[1:datcnt]
    #db.data = data
end



function basgen(datlist::Array, tau, ndim, ires, datcnt, maxbox)
    # Database generator for fet.m function
    # Taehyeun Park; The Cooper Union; EE'15
    
    
    #data = zeros(1,datcnt)
    data = datlist
    trck = 1
    start = 1
    fin = 0
    
    #for ii = 1:length(x)
    #    #if cmp(x[ii], char(32)) || cmp(x[ii], Char(13)) || cmp(x[ii], Char(10)) || cmp(x[ii], Char(26))
    #    if x[ii] == Char(32) || x[ii] == Char(13) || x[ii] == Char(10) || x[ii] == Char(26)
    #        if fin >= start()
    #            data[trck] = parse.(Float64,x[start:fin])
    #            trck = trck + 1
    #            if trck .> 8*floor(datcnt/8)
    #                break
    #            end
    #        end
    #        start = ii + 1
    #    else
    #        fin = ii
    #    end
    #end




    
    delay = collect(0:tau:(ndim-1)*tau)
    
    nxtbox = zeros(maxbox, ndim)
    wherear = zeros(maxbox, ndim)
    datptr = zeros(1,maxbox)
    nxtdat = zeros(1,datcnt)
    
    datmin = minimum(data)
    datmax = maximum(data)
    
    datmin = datmin - 0.01*(datmax - datmin)
    datmax = datmax + 0.01*(datmax - datmin)
    boxlen = (datmax - datmin)/ires
    
    boxcnt = 1
    
    for ii = 1:(datcnt-(ndim-1)*tau)
        target = floor.((data[ii .+ delay] .- datmin)/boxlen)
        runner = 1
        chaser = 0
        
        jj = 1
        while jj <= ndim
            tmp = wherear[Int(runner),jj]-target[jj]
            if tmp < 0
                chaser = runner
                runner = Int(nxtbox[runner,jj])
                if runner != 0
                    continue
                end
            end
            if tmp != 0
               boxcnt = boxcnt + 1
               
               if boxcnt .== maxbox
                   error("Grid overflow, increase number of box count")
               end
               
               for kk = 1:ndim
                   wherear[boxcnt,kk] = wherear[chaser,kk]
               end
               wherear[boxcnt,jj] = target[jj]
               nxtbox[chaser,jj] = boxcnt
               nxtbox[boxcnt,jj] = runner
               runner = Int(boxcnt)
            end
            jj = jj + 1
        end
        nxtdat[ii] = datptr[runner]
        datptr[runner] = ii
    end
    
    used = 0
    for ii = 1:boxcnt
        if datptr[ii] != 0
            used = used + 1
        end
    end
    db = Basegen_db(ndim,ires,tau,datcnt,boxcnt,datmax,datmin,boxlen,datptr[1:boxcnt],nxtbox[1:boxcnt,1:ndim],wherear[1:boxcnt,1:ndim],nxtdat[1:datcnt],data)
    #db.ndim = ndim
    #db.ires = ires
    #db.tau = tau
    #db.datcnt = datcnt
    #db.boxcnt = boxcnt
    #db.datmax = datmax
    #db.datmin = datmin
    #db.boxlen = boxlen
    
    #db.datptr = datptr[1:boxcnt]
    #db.nxtbox = nxtbox[1:boxcnt, 1:ndim]
    #db.wherear = wherear[1:boxcnt, 1:ndim]
    #db.nxtdat = nxtdat[1:datcnt]
    #db.data = data
end


function fet(db, dt, evolve, dismin, dismax, thmax)
    # Computes Lyapunov exponent of given data & parameters; generates output
    # textfile; exact replica of Fortran 77 version of fet()
    # Taehyeun Park; The Cooper Union; EE'15
    
    out = zeros(6)'
    
    ndim = db.ndim
    ires = db.ires
    tau = db.tau
    datcnt = db.datcnt
    datmin = db.datmin
    boxlen = db.boxlen
    
    datptr = db.datptr
    nxtbox = db.nxtbox
    wherear = db.wherear
    nxtdat = db.nxtdat
    data = db.data
    
    delay = 0:tau:(ndim-1)*tau
    datuse = datcnt-(ndim-1)*tau-evolve
    
    its = 0
    SUM = 0
    savmax = dismax
    
    oldpnt = 1
    newpnt = 1
    
    fileID = open("fetout.txt", "w")
    
    goto50 = 1
    while goto50 .== 1
        goto50 = 0
        bstpnt, bstdis, thbest = search(0, ndim, ires, datmin, boxlen, nxtbox, wherear,
            datptr, nxtdat, data, delay, oldpnt, newpnt, datuse, dismin, dismax,
            thmax, evolve)
       
        while bstpnt .== 0
            dismax = dismax * 2
            bstpnt, bstdis, thbest = search(0, ndim, ires, datmin, boxlen, nxtbox, wherear,
                datptr, nxtdat, data, delay, oldpnt, newpnt, datuse, dismin, dismax,
                thmax, evolve)
        end
        
        dismax = savmax
        newpnt = bstpnt
        disold = bstdis
        iang = -1
        
        goto60 = 1
        while goto60 .== 1
            goto60 = 0
            
            oldpnt = oldpnt + evolve
            newpnt = newpnt + evolve
            
            if oldpnt >= datuse
                return out, SUM
            end
            
            if newpnt >= datuse
                oldpnt = oldpnt - evolve
                goto50 = 1
                break
            end
            
            p1 = data[Int.(oldpnt .+ delay)]
            p2 = data[Int.(newpnt .+ delay)]
            disnew = sqrt(sum((p2 - p1).^2))
            
            its = its + 1
    
            SUM = SUM + log(disnew/disold)
            zlyap = SUM/(its*evolve*dt*log(2))
            out = vcat(out, [its*evolve disold disnew zlyap oldpnt-evolve newpnt-evolve])
            #out = [out; its*evolve, disold, disnew, zlyap, (oldpnt-evolve), (newpnt-evolve)];
            
            if iang .== -1
                write(fileID, "$(out[end,1:4]) \n")
                #fprintf(fileID, "#-d\t\t\t#-8.4f\t\t#-8.4f\t\t#-8.4f\n', out[end,1:4]")
            else()
                write(fileID,"$(out[end,1:4]), $iang\n")
                #fprintf(fileID, "#-d\t\t\t#-8.4f\t\t#-8.4f\t\t#-8.4f\t\t#-d\n', [out[end,1:4], iang]")
            end
    
            if disnew <= dismax
                disold = disnew
                iang = -1
                goto60 = 1
                continue
            end
    
            bstpnt, bstdis, thbest = search(1, ndim, ires, datmin, boxlen, nxtbox, wherear, 
                datptr, nxtdat, data, delay, oldpnt, newpnt, datuse, dismin, dismax,
                thmax, evolve)
    
            if bstpnt != 0
                newpnt = bstpnt
                disold = bstdis
                iang = floor(thbest)
                goto60 = 1
                continue
            else()
                goto50 = 1
                break
            end
        end
    end
    close(fileID)
    return out, SUM
end

function solve(prob::LCEProblem, alg::WolfAlgorithm)

    tau = prob.embedded_data.tau
    ndim = prob.embedded_data.dim
    ires = alg.ires
    datcnt = length(prob.timeseries)
    maxbox = alg.maxbox
    datlist = prob.timeseries

    db = basgen(datlist, tau, ndim, ires, datcnt, maxbox)



    dt = alg.dt
    evolve = alg.evolve
    dismin = alg.dismin
    dismax = alg.dismax
    thmax = alg.thmax
    out, SUM = fet(db, dt, evolve, dismin, dismax, thmax)
    lyaps = out[:,4]
    LCEMaxSolution(lyaps[end],lyaps,alg)

end

 
