#include <math.h>
#include <fstream>
#include <iostream>
#include <unordered_set>
#include <sstream>
#include <cstring>
#include <deque>
#include "krfp.h"
#include "minimizers.h"

using namespace std;
typedef uint64_t INT;  

/* Computes the minimizers of a string of length n in O(n) time */
// Implementation storing minimizers in an unordered_set
INT compute_minimizers(  string& text, INT w, INT k, unordered_set<uint64_t> &minimizers )
{
	INT n = text.length();
	INT fp = 0;
	INT smallest_fp = fp;
	
	INT * FP = ( INT * ) malloc( ( n - k + 1  ) *  sizeof( INT ) );
	
	for(INT j = 0; j<k; j++)
		fp =  karp_rabin_hashing::concat( fp, text[j] , 1 );
		
	FP[0] = fp;
	INT pos = 1;
	deque<pair<INT,INT>> min_fp = {};
	
	for(INT j = 1; j<=n-k; j++)
	{
		fp = karp_rabin_hashing::concat( fp, text[j+k-1] , 1);
		fp = karp_rabin_hashing::subtract_k( fp, text[j-1]);
		FP[pos] = fp;
		pos++;
	}	
	
		
	// minimum fp in first window
   	for (INT j = 0; j <= w - k ; j++) 
   	{
 		while ( !min_fp.empty() && FP[j] < min_fp.back().first)
 			min_fp.pop_back();
 			
		min_fp.push_back(std::make_pair(FP[j], j));
    	}
	minimizers.insert( min_fp.at(0).second );
	
	// minium fp in remaining windows
	for( INT i = w-k+1; i<=n-k; i++ )
	{
		while (!min_fp.empty() && min_fp.back().first > FP[i])
			min_fp.pop_back();
	
		min_fp.push_back(std::make_pair(FP[i], i));
		
	
		while( !min_fp.empty() && min_fp.front().second <= i - w + k)
			min_fp.pop_front();
		
		minimizers.insert( min_fp.at(0).second );
	}
	
	free(FP);
	return 0;
}

/* Computes the minimizers of a string of length n in O(n) time */
// Implementation storing minimizers in a vector (with suppressed return)
INT compute_minimizers2(string& text, INT w, INT k)
{
	INT n = text.length();
	INT fp = 0;
	INT smallest_fp = fp;
	vector<INT> minimizers;
	
	INT * FP = (INT *) malloc((n - k + 1) * sizeof(INT));
	
	for(INT j = 0; j<k; j++)
		fp =  karp_rabin_hashing::concat(fp, text[j] , 1);
		
	FP[0] = fp;
	INT pos = 1;
	deque<pair<INT,INT>> min_fp = {};
	
	// find all fps for all k substrings
	for(INT j = 1; j <= n - k; j++)
	{
		fp = karp_rabin_hashing::concat(fp, text[j + k - 1] , 1);
		fp = karp_rabin_hashing::subtract_k(fp, text[j - 1]);
		
		FP[pos] = fp;
		pos++;
	}	
	
	// minimum fp in first window
   	for (INT j = 0; j <= w - k ; j++) 
   	{
 		while ( !min_fp.empty() && FP[j] < min_fp.back().first)
 			min_fp.pop_back();
				
		min_fp.push_back(std::make_pair(FP[j], j));
    	}
	minimizers.push_back( min_fp.at(0).second );
	
	// minium fp in remaining windows
	for( INT i = w - k + 1; i <= n - k; i++)
	{
		while (!min_fp.empty() && min_fp.back().first > FP[i])
			min_fp.pop_back();

		min_fp.push_back(std::make_pair(FP[i], i));
		
	
		while(!min_fp.empty() && min_fp.front().second <= i - w + k)
		{
			min_fp.pop_front();
		}	
		if(minimizers.back() != min_fp.at(0).second)
			minimizers.push_back(min_fp.at(0).second);
	}
	
	free(FP);
	return 0;
}
