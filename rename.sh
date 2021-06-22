a=1
for i in *Pass.withcoords; do
  new=$(printf "%02d.Pass.withcoords" "$a") #02 pad to length of 2
  cp -i -- "$i" "$new"
  let a=a+1
done
