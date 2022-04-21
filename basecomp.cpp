#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <vector>
#include <string>


std::vector<std::vector<unsigned int>> compute(unsigned int refLen, std::vector<std::string> reads, std::vector<unsigned int> starts, std::vector<std::vector<std::pair<unsigned int, unsigned int>>> ctuples)
{
    // 2D array containing counts of A, C, G, T, N (Illumina unknown base) and DS (deletions+skips) present in the reads, for each position in the reference
    std::vector<std::vector<unsigned int>> baseCounts(refLen, std::vector<unsigned int>(6));

    // Iterate through the reads
    for (unsigned int i = 0; i < reads.size(); i++)
    {
        // The aligned portion of the read (soft clipped bases are excluded)
        auto &read = reads[i];

        // Sequence of operations describing how the read is aligned to the reference
        // Each pair is of the form (operation, length of operation)
        auto &tups = ctuples[i];
        
        // Index of where the read starts in the reference
        unsigned int refPos = starts[i];

        // Index for position in (the aligned portion of) the current read
        unsigned int seqPos = 0;

        for (auto &tup : tups)
        {
            // Integer designating the operation
            auto &operation = tup.first;

            // Length of operation
            auto &opLen = tup.second;

            // 0 : Matched (equal or not)
            // 7 : All bases equal to ref
            // 8 : All bases not equal to ref
            if ((operation == 0) || (operation == 7) || (operation == 8))
            {
                // Iterate through the read for the given length of the match
                for (unsigned int j = 0; j < opLen; j++)
                {
                    switch(read[seqPos + j])
                    {
                        case 'A' : baseCounts.at(refPos + j).at(0) += 1; break;
                        case 'C' : baseCounts.at(refPos + j).at(1) += 1; break;
                        case 'G' : baseCounts.at(refPos + j).at(2) += 1; break;
                        case 'T' : baseCounts.at(refPos + j).at(3) += 1; break;
                        case 'N' : baseCounts.at(refPos + j).at(4) += 1; break;
                    }
                }

                // Update the current positions in the reference and the read
                refPos += opLen;
                seqPos += opLen;
            }

            // 1 : Insertion in read
            // Because it is an insertion, move through the insert without moving further in the reference
            else if (operation == 1)
                seqPos += opLen;

            // 2 : Deletion in read
            // 3 : Skip in read
            // Because it is a deletion/skip, move through the reference without moving further in the read
            else if ((operation == 2) || (operation == 3))
            {
                // Count the positions where deletions/skips occur
                for (unsigned int j = 0; j < opLen; j++)
                    baseCounts.at(refPos + j).at(5) += 1;

                refPos += opLen;
            }

            // Operations 4, 5 and 6 are soft clipping, hard clipping and padding respectively
            // None of these affect the read.query_alignment_sequence                    
            // Operation 9 is the 'back' operation which seems to be basically unheard of
            // TODO: can it be ignored?
        }
    }
    return baseCounts;
}


PYBIND11_MODULE(basecomp, m)
{
    m.def("compute", &compute);
}
