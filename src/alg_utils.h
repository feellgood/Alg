#ifndef ALG_UTILS_H
#define ALG_UTILS_H

#include<type_traits>
#include <vector>
#include <iostream>
#include <algorithm>

#include "alg_coeff.h"

namespace alg
{

/** overloaded operator for vector printing */
inline std::ostream & operator<<(std::ostream & flux, std::vector<double> const& v) 
{
std::for_each(v.begin(),v.end(), [&flux](const double& x) { flux << x << " "; });
return flux;
}

enum precondType {NONE, DIAG};

template <typename T>
class CSR_mat
	{
	public:
	/**
	constructor : memory allocation only
	*/
		CSR_mat(const int _N, const int nb_coeff) :N(_N) 
			{
			I = new int[N+1];
			I[N] = nb_coeff;
			J = new int[nb_coeff];
			val = new T[nb_coeff];
			}
		
		/**
constructor : build the CSR sparse matrix from the vector of m_coeff C ordered
dim_C is the dimension of the ending space of the application C = number of lines.
Zero index based.
*/
		CSR_mat(const std::vector<alg::m_coeff> C, const size_t dim_C) : N(dim_C)
		{
		size_t nb_coeff = C.size();
		I = new int[N+1];
		I[N] = nb_coeff;
		J = new int[nb_coeff];
		val = new T[nb_coeff];

		size_t i(0);
		for(size_t k=0;k<nb_coeff;k++)
			{
			val[k] = (T) C[k].getVal();
			J[k]=C[k]._j;
			if((C[k]._i == i )&&(i < dim_C)) 
				{
				I[i] = k;  
				i++;
				}
			}
		}
		
		/**
		constructor copy : we do not need it, so we define it deleted.
		*/
		CSR_mat(const CSR_mat&) = delete;
		
		/**
		assignment operator= :  we do not need it, so we define it deleted.
		*/
		CSR_mat& operator=(const CSR_mat&) = delete;
		
		~CSR_mat()
			{
			delete[] I;
			delete[] J;
			delete[] val;
			}
		
		int *I;
		int *J;
		T *val;
		const int N;// number of line
	};

/** return square(x) */
template <typename T>
inline T sq(T x) {return x*x;}

/** multiplication of a sparse matrix with a dense vector ; CSR sparse matrix ; return : y , y must be zero initialized  */
template <typename T>
void multCSR_MatVect(CSR_mat<T> const& A,T *x,T *y)
{
for(int i=0;i<A.N;i++)
	for(int k=A.I[i];k<A.I[i+1];k++)
		{ y[i] += A.val[ k ]*x[A.J[k]]; }
}

/**
build diagonal preconditionner of CSR matrix ; all diagonal coefficients must be non zero. If a zero is on the diagonal diag is left unmodified.
diagDim might be different from the dimension of the matrix
*/
template <typename T>
void build_diagPrecond_CSRmat(CSR_mat<T> const& A, T * diag,const int diagDim)
{
for(int i=0;i<diagDim;i++)
	for(int k=A.I[i];k<A.I[i+1];k++)
		{ 
		if ((A.J[k] == i)&&(A.val[k]!=0.0))
			{ diag[i] = 1.0/A.val[k]; }
		}
}


/** 
in place left multiplication of a CSR sparse matrix by a diagonal matrix stored in diag
 */
template <typename T>
void leftMult_diagPrecond_CSRmat(CSR_mat<T> & A, T * diag,const int diagDim)
{
// diagDim might be inferior to N
if(diagDim>A.N)
	{ exit(1); }
else
	{
	for(int i=0;i<diagDim;i++)
		for(int k=A.I[i];k<A.I[i+1];k++)
			{ A.val[k] *= diag[i]; }
	}
}

/**
computes norm2(A x - rhs) 
*/
template <typename T>
T check_sol(alg::CSR_mat<T> const& A,T *x,T *rhs)
{
T result(0);
const int DIM_y = A.N;
T y[DIM_y];

for(int k =0;k<DIM_y;k++) {y[k]= 0;} //initialization to zero

alg::multCSR_MatVect<T>(A,x,y);

for(int k =0;k<DIM_y;k++) 
	{ result += alg::sq<T>( y[k] - rhs[k]); }

return sqrt(result);
}

/**
node template for btree
*/
template<typename T> 
struct node_t{
    T *data;
	node_t *left;
	node_t *right;
};

/**
btree for write sparse matrix
*/
template<typename T> 
class btree{
public:
	/** constructor */
	btree(): root(NULL) {}
	
	/** destructor */
	~btree() { destroy_tree(root); }


	/**
	constructor copy : we do not need it, so we define it deleted.
	*/
	btree(const btree&) = delete;
		
	/**
	assignment operator= :  we do not need it, so we define it deleted.
	*/
	btree& operator=(const btree&) = delete;

	/**
	data inserter
	*/
	void insert(T data);
	
	node_t<T> *search(size_t _i)
		{ return search(_i, root); }
	
	/** printing function */
	void inorder_print()
		{ inorder_print(root); std::cout << "\n"; }
	
	/** printing function */
	void postorder_print()
		{ postorder_print(root); std::cout << "\n"; }
	
	/** printing function */
	void preorder_print()
		{ preorder_print(root); std::cout << "\n"; }
    
    	/** std::vector inserter  */
    	void inorder_insert(std::vector<T> &v)
    		{ inorder_insert(root, v); }

private:
	void destroy_tree(node_t<T> *leaf);
	void insert(T data, node_t<T> *leaf);
	node_t<T> *search(size_t _i, node_t<T> *leaf);
	void inorder_print(node_t<T> *leaf);
	void postorder_print(node_t<T> *leaf);
	void preorder_print(node_t<T> *leaf);
    	void inorder_insert(node_t<T> *leaf, std::vector<T> &v);
	
	node_t<T> *root;
};


/** 
public data inserter
*/
template<typename T> 
void btree<T>::insert(T data)
{
	if(root != NULL)
		{ insert(data, root); }
	else
		{
		root = new node_t<T>;
        	if (!root) exit(1);
        	root->data = new T(data._i, data.getVal());
		root->left = NULL;
		root->right = NULL;
		}
}


/**
private function to destroy the tree and dealloc mem, recursive
*/
template<typename T> 
void btree<T>::destroy_tree(node_t<T> *leaf){
	if(leaf != NULL){
		destroy_tree(leaf->left);
		destroy_tree(leaf->right);
        	delete leaf->data;
		delete leaf;
	}
}


/**
private inserter
*/
template<typename T> 
void btree<T>::insert(T data, node_t<T> *leaf){

	if (data._i == leaf->data->_i){
       leaf->data->add(data.getVal());
       return;
	   }

	if (data._i < leaf->data->_i){
	   if (leaf->left != NULL){
		  insert(data, leaf->left);
		}
       else{
		  leaf->left = new node_t<T>; 
          if (!leaf->left) exit(1);
          T* data_ptr = new T(data._i, data.getVal()); 
          if (!data_ptr) exit(1);
          leaf->left->data  = data_ptr;
		  leaf->left->left  = NULL;
		  leaf->left->right = NULL;
		  }
       return;
	   }

    if (data._i > leaf->data->_i){
	   if (leaf->right != NULL){
		  insert(data, leaf->right);
		  } 
       else{
		  leaf->right = new node_t<T>;
          if (!leaf->right) exit(1);
          T* data_ptr = new T(data._i, data.getVal()); 
          if (!data_ptr) exit(1);
          leaf->right->data  = data_ptr;
		  leaf->right->right = NULL;
		  leaf->right->left  = NULL;
		  }
       return;
	   }
}


/** 
private leaf search, recursive
*/
template<typename T> 
node_t<T> *btree<T>::search(size_t _i, node_t<T> *leaf){
	if(leaf != NULL){
		if(_i == leaf->data->_i){
			return leaf;
		}
		if(_i < leaf->data->_i){
			return search(_i, leaf->left);
		}else{
			return search(_i, leaf->right);
		}
	}else{
		return NULL;
	}
}

/** private printing function */
template<typename T>
void btree<T>::inorder_print(node_t<T> *leaf){
	if(leaf != NULL){
		inorder_print(leaf->left);
		std::cout << "{" << leaf->data->_i<< ":"<<leaf->data->getVal() << "},";
		inorder_print(leaf->right);
	}
}

/** private printing function */
template<typename T>
void btree<T>::postorder_print(node_t<T> *leaf){
	if(leaf != NULL){
		inorder_print(leaf->left);
		inorder_print(leaf->right);
		std::cout << "{" << leaf->data->_i<< ":"<<leaf->data->getVal() << "},";
	}
}

/** private printing function */
template<typename T>
void btree<T>::preorder_print(node_t<T> *leaf){
	if(leaf != NULL){
		std::cout << "{" << leaf->data->_i<< ":"<< leaf->data->getVal() << "},";
		inorder_print(leaf->left);
		inorder_print(leaf->right);
	}
}

/** private inserter for vector */
template<typename T> 
void btree<T>::inorder_insert(node_t<T> *leaf, std::vector<T> &v){
	if(leaf != NULL)
		{
		inorder_insert(leaf->left, v);
        	alg::v_coeff data{leaf->data->_i, leaf->data->getVal()};
        	v.push_back(data);
		inorder_insert(leaf->right, v);
		}
}

} // end namespace alg

#endif

