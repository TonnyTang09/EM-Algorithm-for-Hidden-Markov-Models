# EM-Algorithm-for-Hidden-Markov-Models
The special feature of this program is the construction of the parameter estimation method, i.e., 
the EM algorithm, for Hidden Markov Models whose output distributions are different
P-dimensional normal distributions in each state.Since the EM algorithm is an initial value sensitive iterative algorithm, 
this program first clusters the data to get the initial parameters, 
and then uses the EM algorithm for more detailed learning. 
In the construction process, this program utilizes the idea of LogSumExp, 
which successfully avoids the problem of data underflow. 
This program also implements the Viterbi algorithm to make predictions about the hidden Markov chain states. 
Finally the program iteration speed is accelerated according to the parameter expansion method proposed by Rubin in his 1998 article.

