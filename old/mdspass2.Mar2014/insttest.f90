
  double precision ar(3,3),e(3),ei(3),evi(3,3),evr(3,3),work(12)
  nn=3
  ar(1,1)=0.0d0; ar(1,2)=1.0; ar(1,3)=1.0;
  ar(2,1)=1.0d0; ar(2,2)=0.0; ar(2,3)=1.0;
  ar(3,1)=1.0d0; ar(3,2)=1.0; ar(3,3)=0.0;
  

  call dgeev('N','V',nn,ar,nn,e,ei,evi,nn,evr,nn,work,nn*4,INFO)

  print *,evr(1,3)
  print *,evr(2,3)
  print *,evr(3,3)

  end
