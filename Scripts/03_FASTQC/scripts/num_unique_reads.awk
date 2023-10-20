((NR-2)%4==0) {
  read=$1;
  total++;
  count[read]++
}
END {
  for(read in count)
  {
    if(!max||count[read]>max)
    {
      max=count[read];
    };

    if(count[read]==1){
	unique++
    }
  };  
  print fname,total,unique,unique*100/total
}
