#include<math.h> // sqrt
#include<vector>
#include<stdio.h>

using namespace std;

struct I{
	// imaginary numbers
	float real,imag;
	I(float r,float i=0):real(r),imag(i){}
	I operator+(I a){return (I){real+a.real,imag+a.imag};}
	I operator*(I a){return (I){real*a.real-imag*a.imag,real*a.imag+a.real+imag};}
	I operator+=(I a){return*this=*this+a;}
	float sqabs(){return real*real+imag*imag;}
	float abs(){return sqrt(sqabs());}
};


I sqrt(I z){
	float x=z.real,y=z.imag;

	if(x == 0){
		float t = sqrt(abs(y)/2);
		return (I){t,y<0?-t:t};
	}else{
		float t = sqrt(2*(z.abs()+abs(x)));
		float u = t / 2;
		return x > 0
		? (I){u, y/t}
		: (I){abs(y)/t,y<0?-u:u};
	}
}


struct Filter{
	int n;
	I*matrix;
	Filter(int n):n(n){matrix=(I*)malloc(n*n*4*sizeof(I));} // 2^n^2 lol
	//Filter(){} // XXX initializer list
};

void matmult(int s,I*in,I*mat){
	I*res=(I*)malloc(s*sizeof(I));
	for(int x=0;x<s;x++)
	for(int y=0;y<s;y++)
	res[x]+=in[y]*mat[x*s+y];

	for(int i=0;i<s;i++)
	in[i]=res[i];

	// returns in `in'
}

class Qb{public:
	int n; // how many there are
	I*data;

	Qb(int size,char init=0):n(size){
		data=(I*)malloc((1<<n)*sizeof(I));
		for(int i=0;i<(1<<n);i++)
			data[i]=(I){init-1?0:(float)i,0};

		if(init==0)
			data[(1<<n)-1]=(I){1,0};
	}

	void set(int index,bool state){
		// forall
		//  if it represents the state we want
		//   figure out the complement state
		//   the probability of this state is the euclidean of them both
		//   the probability of the complement is zero
		for(int i=0;i<(1<<n);i++){
			if(i&&(1<<index) == state){
				int c=state?(i&!(1<<index)):(i|(1<<index));
				data[i]=sqrt(data[i]*data[i]+data[c]*data[c]); // XXX idek, should be right
				data[c]=(I){0,0};
			}
		}
	}

	float observe(int index){ // XXX
		// iterate thru all possibilities
		//  figure out whether they increment or decrement probability
		//  increment/decrement probability accordingly
		// either pick random or most likely
		// set()
		// return

		float prob=0;
		for(int i=0;i<(1<<n);i++)
			if(i&(1<<index)){
				printf("p+%.3f @ %i\n",data[i],i);
				prob+=data[i]*data[i]; // XXX i am not sure what it needs here
			}

		return prob;

		// actually, given the cached data, it might be easier to set it ourselves
		return 1;
	}

	void apply(Filter f,vector<int>indices){
		// we iterate over all possible arrangements of other qubits

		// for 0<i<n-f.n
		//  shift bits around to get all but the indexed XXX fuck this shit hard
		//  arrange indices in a new vector
		//  arrange input in a new vector
		//  multiply vector w/ filter
		//  arrange output back into data (according to index vector)

		// XXX ok uh, this might be a bit more complicated than a multiplication
		// . . . or is it

		printf("%i / %i qubits\n",f.n,n);

		for(int i=0;i<(1<<n-f.n);i++){
			//printf("combination %x\n",i);

			int mask = i;
			for(int a : indices){
				/*
				printf("m: %i %i %i %i\n",
					mask,a,
					mask>>a<<a+1,
					(mask&((1<<a)-1))
				);
				*/
				//mask = (mask<<1)&!((1<<a)-1) | (mask&((1<<a)-1));
				mask = (mask>>a<<a+1) | (mask&((1<<a)-1));
			}

			printf("mask %x\n",mask);

			vector<int>ind,tmp;

			ind.push_back(mask);

			for(int a : indices){
			for(int b : ind)
				tmp.push_back(b|(1<<a));
			for(int b : tmp)
				ind.push_back(b);
			tmp.clear();
			}

			printf("indices: ");
			for(int a:ind)printf("%i ",a);
			putchar(10);

			I*v=(I*)malloc((1<<f.n)*sizeof(I));
			for(int i=0;i<(1<<f.n);i++)
				v[i]=data[ind[i]];
			matmult(1<<f.n,v,f.matrix);
			for(int i=0;i<(1<<f.n);i++)
				data[ind[i]]=v[i];
		}
	}

	void print(){
		for(int i=0;i<(1<<n);i++){
			printf("probability for ");
			//for(int j=0;j<n;j++)
			for(int j=n-1;j!=-1;j--)
				putchar(i&(1<<j)?'1':'0');
			printf(" : %g\n",i,data[i]);
		}
		putchar(10);
	}
};



int main(){
	Qb comp=Qb(3,0);
	comp.print();

	//comp.data[2]=comp.data[6]=0.7071067811865476;comp.data[7]=0;

	// CNOT
	Filter cnot=Filter(2);
	cnot.matrix[0]=I(1);
	cnot.matrix[5]=I(1);
	cnot.matrix[11]=I(1);
	cnot.matrix[14]=I(1);

	// Hadamard
	float insqrt=1/sqrt(2);
	Filter had=Filter(1);
	had.matrix[0]=I( insqrt);
	had.matrix[1]=I( insqrt);
	had.matrix[2]=I( insqrt);
	had.matrix[3]=I(-insqrt);

	printf("HADAMARD(2)\n");
	comp.apply(had,{2});
	comp.print();

	printf("CNOT(0,1)\n");
	comp.apply(cnot,{0,1});
	comp.print();

	return 0;
}

/*
to test:
- observation
- initialization
- basic gates and application

to define:
- observation
- initialization
- basic gates

gates:
at the very least, hadamard/1 and controlled-v/2, used in qbf. though, id prefer to have more
hadamard: [1,1,1,-1]/sqrt(2)
controlled-v: multiplies by i if both 1.
swap/2 can be useful: [1000 0010 0100 0001]
toffoli/3 (or ccnot) is a not gate w/ two controls
*/
