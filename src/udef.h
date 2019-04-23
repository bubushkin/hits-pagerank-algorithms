/*
 * udef.h
 *
 *  Created on: Apr 22, 2019
 *      Author: iskandar
 */

#ifndef UDEF_H_
#define UDEF_H_

#define CMD_ALGO argv[1]
#define CMD_ITERATIONS argv[2]
#define CMD_INITVAL argv[3]
#define CMD_INPUTFILE argv[4]

#define NEWLINE std::cout << std::endl;

#define HELP()  cout << "HELP USAGE: "<< endl; \
                cout << argv[0x0] << " ALGO_OPTION ITERATIONS INITVALUE INPUT_FILE"<< endl; \
                cout << "OPTIONS:" << endl; \
                cout << "\t-h\t:run HITS Algorithm" << endl; \
                cout << "\t-p\t:run Pagerank Algorithm" << endl; \


#endif /* UDEF_H_ */
