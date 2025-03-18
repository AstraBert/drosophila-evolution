echo "will cite" | parallel --citation >/dev/null 2>&1
parallel --bar -j 100 bash ::: ../../shell/wget_data/*.sh
