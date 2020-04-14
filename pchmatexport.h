/*
________________________

Author: Deniz Bilgili
Date  : April 12th, 2020
________________________

MIT License

Copyright (c) 2020 Deniz Bilgili

Permission is hereby granted, free of charge, to any person obtaining a copy of 
this software and associated documentation files (the "Software"), to deal in the 
Software without restriction, including without limitation the rights to use, 
copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the 
Software, and to permit persons to whom the Software is furnished to do so, 
subject to the following conditions:

The above copyright notice and this permission notice shall be included in all 
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR 
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS 
FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR 
COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN 
AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION 
WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

#ifndef PCHMATEXPORT_H
#define PCHMATEXPORT_H

#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <map>
#include <Eigen/Dense>

namespace PCH {

enum class MATRIX : bool {
	MASS = true,
	STIF = false
};

size_t dofs_before(size_t nodeID,
	const std::map<size_t, size_t>& dofmap) {
	/*
	Given the ID of a node, nodeID, returns the sum of the number of DOFs
	belonging to the nodes before the given node. dofmap is a hash table that maps
	node IDs to the total number of DOFs assigned to that node. This function goes
	through dofmap and sums the total number of DOFs for each node before the
	given node. This is later used together with build_dofmap() for easy indexing
	of the system matrix elements during the extraction of matrices from PCH file.
	*/
	std::map<size_t, size_t>::const_iterator itr;

	size_t dofs = 0;

	for (itr = dofmap.cbegin(); itr != dofmap.cend(); itr++) {
		if (itr->first == nodeID) break;
		dofs += itr->second;
	}

	return dofs;
}

std::map<size_t, size_t> build_dofmap(std::string pchAddress) {
	/*
	Returns a hash table that maps the IDs of the nodes to the total number of
	DOFs the nodes have. For example, if the node with the ID 82 has 3 DOFs
	defined, then its entry on the table would be {82 : 3}. Dofmap is used for
	easy indexing of the system matrix elements during the extraction of matrices
	from PCH files.
	*/
	std::ifstream ifs(pchAddress);
	std::string lineRaw;

	size_t numDofs; // compare this later with the sum of values in dofmap
	std::map<size_t, size_t> dofmap;

	while (std::getline(ifs, lineRaw)) {
		std::istringstream iss(lineRaw);
		std::string lineToken;

		iss >> lineToken;

		if (lineToken == "SPOINT") {
			std::string modeIDStr;
			iss >> modeIDStr;
			size_t modeID = std::stoi(modeIDStr);

			dofmap.insert(std::pair<size_t, size_t>(modeID, 1));

			continue;
		}

		if (lineToken == "DMIG") {
			std::vector<std::string> dmigInfo;

			while (iss >> lineToken) {
				dmigInfo.push_back(lineToken);
			}

			numDofs = std::stoi(dmigInfo.back());

			break;
		}
	}

	while (std::getline(ifs, lineRaw)) {
		std::istringstream iss(lineRaw);
		std::string lineToken;

		iss >> lineToken;

		if (lineToken == "DMIG") break;
		if (lineToken != "DMIG*") continue;

		std::string matrixType;
		std::string nodeIDStr;
		std::string nodeDofStr;

		iss >> matrixType;
		iss >> nodeIDStr;
		iss >> nodeDofStr;

		if (nodeDofStr == "0") continue;

		size_t nodeID = std::stoi(nodeIDStr);
		size_t nodeDof = std::stoi(nodeDofStr);

		dofmap[nodeID] = nodeDof;
	}

	return dofmap;
}

Eigen::MatrixXd read_matrix(std::string pchAddress, MATRIX matrixType) {
	/*
	Extracts the requested system matrix (mass or stiffness) from the given PCH
	file and reads it into the memory. Uses build_dofmap() to generate a Dofmap.
	
	The total number of DOFs is read from the DMIG line and a DOFxDOF square
	matrix is zero initialized. Below is an example how the zero initialized
	system matrix is filled:

	DMIG*           KAAX                 546       3
	*                    547       2         5.26D3

	The row of the matrix element 5.26D3 is equal to the number of DOFs before
	the node 546, plus the current DOF of the node 546, which is 3 in this
	example. So the row is found as below:
					dofsBefore = dofs_before(546, dofmap);
					globalRow = dofsBefore + 3 - 1;
	Similarly, the column of the matrix element 5.26D3 is equal to the number of
	DOFs before the node 547, plus the current DOF of the node 547, which is 2 in
	this example. So the column is found as below:
					dofsBefore = dofs_before(547, dofmap);
					globalCol = dofsBefore + 2 - 1;
	Finally, the corresponding element of the system matrix is updated as below:
					systemMatrix(globalRow, globalCol) = 5.26D3
	Since the system matrix is symmetric, we also do below:
					systemMatrix(globalCol, globalRow) = 5.26D3
	*/
	std::string userMatrixType = "KAAX";
	if (matrixType == MATRIX::MASS) userMatrixType = "MAAX";

	auto dofmap = build_dofmap(pchAddress);
	auto numDofs =
		dofs_before(dofmap.rbegin()->first, dofmap) + dofmap.rbegin()->second;

	Eigen::MatrixXd sysMatrix(numDofs, numDofs);
	sysMatrix.setZero();

	std::ifstream ifs(pchAddress);
	std::string lineRaw;

	std::string fileMatrixType;
	std::string nodeIDStr;
	std::string nodeDofStr;

	size_t nodeID;
	size_t nodeDof;

	size_t globalRow = 0;
	size_t globalCol = 0;

	bool doProcess = false;

	while (true) {
		std::getline(ifs, lineRaw);
		if (ifs.eof()) break;

		std::istringstream iss(lineRaw);
		std::string lineToken;

		iss >> lineToken;

		if (lineToken == "DMIG*") {
			doProcess = true;

			iss >> fileMatrixType;
			iss >> nodeIDStr;
			iss >> nodeDofStr;

			if (fileMatrixType != userMatrixType) {
				doProcess = false;
				continue;
			}

			nodeID = std::stoi(nodeIDStr);
			nodeDof = std::stoi(nodeDofStr);

			auto dofsBefore = dofs_before(nodeID, dofmap);
			globalRow = dofsBefore + nodeDof - 1;
		}

		if (lineToken == "*" && doProcess) {
			std::string matrixElemStr;

			iss >> nodeIDStr;
			iss >> nodeDofStr;
			iss >> matrixElemStr;

			for (auto& c : matrixElemStr) {
				if (c == 'D') c = 'E';
			}

			nodeID = std::stoi(nodeIDStr);
			nodeDof = std::stoi(nodeDofStr);
			double matrixElem = std::stod(matrixElemStr);

			auto dofsBefore = dofs_before(nodeID, dofmap);
			globalCol = dofsBefore + nodeDof - 1;

			sysMatrix(globalRow, globalCol) = matrixElem;
			sysMatrix(globalCol, globalRow) = matrixElem;
		}
	}

	return sysMatrix;
}

}

#endif // !PCHMATEXPORT_H
