find . -mindepth 3 -depth -type d -name 'gen*' -exec rm -rf '{}' \;
find . -name 'summary*' -delete
for file in readme.err readme.log run_lamp_*.err run_lamp_*.log; do
   if [ -e $file ]; then rm $file; fi
done
