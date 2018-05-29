IC=3            #{Entropy Gaussian, Entropy Shock, Sod Shocktube}
Amp=3
T=0.1
L=1.5

CFLs=( 0.45 )
Ns=( 128 256 512 1024 2048 4000 8000 16000 32000 64000 )
#Ns=( 51 )

Beta=0.9                        #for artifical viscosity
out_file=200

enable_adjoint='True'           #{True,False}
enable_viscosity='True'         #{True,False}
enable_f_viscosity='True'      #{True,False}
#enable_diagonalization='Split'  #{True,False,Split} 
enable_acc='False'

ForTime=Euler       #{Euler,RK3}
BackTime=Euler      #{Euler,RK3}
forflux=GLF         #{GLF,AV,LLF}
backflux=GLF        #{GLF,AV}
SI=FV               #{WENO,FV}

#ed=( 'True' 'False' 'Split' )
ed=( 'True' 'False' )
#ed=( 'True' )
rm -rf ./out/* ./O/*

for enable_diagonalization in "${ed[@]}"; do
for N in "${Ns[@]}"; do
for CFL in "${CFLs[@]}"; do
    case='weno-'$T'T-'$N'N-'$CFL'CFL'
    echo Case $case
    cat ./Input/weno.temp   | sed -e s/XXACC/$enable_acc/ -e s/XXFVISC/$enable_f_viscosity/ -e s/XXDIAG/$enable_diagonalization/ -e s/XXSI/$SI/ -e s/XXOUT/$out_file/g -e s/XXIC/$IC/g -e s/XXN/$N/ -e s/XXBETA/$Beta/ -e s/XXCFL/$CFL/ -e s/XXT/$T/ -e s/XXAMP/$Amp/ -e s/XXL/$L/ -e s/XXFT/$ForTime/ -e s/XXBT/$BackTime/ -e s/XXFFLUX/$forflux/ -e s/XXAFLUX/$backflux/ -e s/XXADJ/$enable_adjoint/ -e s/XXVISC/$enable_viscosity/ > ./Input/weno.in
    rm -rf ./D/* 
    #./weno3.x  > './out/'$case'_3.out'
    ./weno5.x  > './out/'$case'_5.out'
done
done

rm -rf ./O_$enable_diagonalization
cp -r ./O ./O_$enable_diagonalization
done
