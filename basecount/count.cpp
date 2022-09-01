#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <vector>
#include <string>


std::vector<std::vector<unsigned int>> bcount(
    const unsigned int refLen,
    const unsigned int minBaseQuality,
    const std::vector<std::string> reads,
    const std::vector<std::vector<unsigned int>> qualities,
    const std::vector<unsigned int> starts, 
    const std::vector<std::vector<std::pair<unsigned int, unsigned int>>> ctuples
    )
{
    // 2D array that, for each position in the reference, counts A, C, G, T, deletions+skips and N (Illumina unclear basecall) present in the reads
    std::vector<std::vector<unsigned int>> baseCounts(refLen, std::vector<unsigned int>(6));

    const unsigned int numReads = reads.size();

    // Iterate through the reads, starts, and ctuples
    for (unsigned int i = 0; i < numReads; i++)
    {
        // Read with soft clipped bases excluded
        const std::string &read = reads[i];

        // Qualities for the read, with soft clipped bases excluded
        const std::vector<unsigned int> &quals = qualities[i];

        // Sequence of operations describing how the read is aligned to the reference
        // Each pair is of the form (operation, length of operation)
        const std::vector<std::pair<unsigned int, unsigned int>> &tups = ctuples[i];
        
        // Current position in the reference
        unsigned int refPos = starts[i];

        // Current position in the read
        unsigned int readPos = 0;

        for (const auto &tup : tups)
        {
            // Integer designating the operation
            const auto &operation = tup.first;

            // Length of operation
            const auto &opLen = tup.second;

            // 0 : Ref and read are matched (bases equal or not)
            // 7 : Ref and read have equal bases
            // 8 : Ref and read have different bases
            if ((operation == 0) || (operation == 7) || (operation == 8))
            {
                // Iterate through the read for the given length of the operation
                for (unsigned int j = 0; j < opLen; j++)
                {
                    if (quals[readPos] >= minBaseQuality)
                    {
                        switch(read[readPos])
                        {
                            case 'A' : baseCounts.at(refPos).at(0) += 1; break;
                            case 'C' : baseCounts.at(refPos).at(1) += 1; break;
                            case 'G' : baseCounts.at(refPos).at(2) += 1; break;
                            case 'T' : baseCounts.at(refPos).at(3) += 1; break;
                            case 'N' : baseCounts.at(refPos).at(5) += 1; break;
                        }
                    }
                    readPos += 1;
                    refPos += 1;
                }
            }

            // 1 : Insertion in read
            // Move through the inserted portion of read, without moving in the reference
            else if (operation == 1)
                readPos += opLen;

            // 2 : Deletion in read
            // 3 : Skip in read
            // Because it is a deletion/skip, move through the reference without moving further in the read
            else if ((operation == 2) || (operation == 3))
            {
                // Count the positions where deletions/skips occur
                for (unsigned int j = 0; j < opLen; j++)
                {
                    baseCounts.at(refPos).at(4) += 1;
                    refPos += 1;
                }
                // TODO: Should deletions be exempt from quality scoring?
                // They are assigned a score in samtools pileup but it is just the nearest quality of a non-deleted base
            }

            // Operations 4, 5 and 6 are soft clipping, hard clipping and padding respectively
            // None of these affect the read.query_alignment_sequence passed in from pysam                    
            // Operation 9 is the 'back' operation which seems to be basically unheard of
            // TODO: can it be ignored?
        }
    }
    return baseCounts;
}


PYBIND11_MODULE(count, m)
{
    m.def("bcount", &bcount);
}
