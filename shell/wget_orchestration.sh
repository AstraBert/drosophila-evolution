cnt=0
mkdir wget_data/

while IFS= read -r url; do
    ((cnt++))
    echo "wget "$url"" > wget_data/${cnt}.sh
done < ../data/bam_download_link.txt
