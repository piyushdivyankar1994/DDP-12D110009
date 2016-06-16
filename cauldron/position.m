A=[0 0 0];

for i1=0:4
	temp1=A+[i1 0 0];
	for i2=0:4
		temp2=temp1+[0 i2 0];
		for i3=0:4
			temp3=temp2+[0 0 i3];
			dist1(25*i1+5*i2+i3+1,1:3)=temp3;
			dist1(25*i1+5*i2+i3+1,4)=norm(temp3);
		endfor
	endfor
endfor		

A=[0.5 0.5 0];

for i1=0:4
	temp1=A+[i1 0 0];
	for i2=0:4
		temp2=temp1+[0 i2 0];
		for i3=0:4
			temp3=temp2+[0 0 i3];
			dist2(25*i1+5*i2+i3+1,1:3)=temp3;
			dist2(25*i1+5*i2+i3+1,4)=norm(temp3);
		endfor
	endfor
endfor		

A=[0 0.5 0.5];

for i1=0:4
	temp1=A+[i1 0 0];
	for i2=0:4
		temp2=temp1+[0 i2 0];
		for i3=0:4
			temp3=temp2+[0 0 i3];
			dist3(25*i1+5*i2+i3+1,1:3)=temp3;
			dist3(25*i1+5*i2+i3+1,4)=norm(temp3);
		endfor
	endfor
endfor		

A=[0.5 0 0.5];

for i1=0:4
	temp1=A+[i1 0 0];
	for i2=0:4
		temp2=temp1+[0 i2 0];
		for i3=0:4
			temp3=temp2+[0 0 i3];
			dist4(25*i1+5*i2+i3+1,1:3)=temp3;
			dist4(25*i1+5*i2+i3+1,4)=norm(temp3);
		endfor
	endfor
endfor	

sorted_dist=vertcat(dist1,dist2,dist3,dist4);
sorted_dist=sortrows(sorted_dist,4);

