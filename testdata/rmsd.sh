
file=$1

theoretical_col=${2:-"3"}
exp_col=${3:-"4"}

grep -v "\(^\s*#\|^\s*$\)" $file | awk '

    BEGIN{
        sum_sd = 0
    }

    {
        sum_sd += ($'$theoretical_col'-$'$exp_col')*($'$theoretical_col'-$'$exp_col')
    }
    
    END{
        rmsd = sqrt(sum_sd / NR)

        printf "\n\t% lf\n\n", rmsd
           
    }'
