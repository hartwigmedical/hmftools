package com.hartwig.hmftools.viridian.reference;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import com.hartwig.hmftools.common.utils.file.DelimFileReader;
import com.hartwig.hmftools.viridian.common.UserInputError;

// Per-virus-genome reference data.
public record VirusInfo(
        String contigName,
        String virusName,
        OncologyGroup oncologyGroup
)
{
    public static List<VirusInfo> load(String tsvFile)
    {
        List<VirusInfo> result = new ArrayList<>();
        Set<String> contigs = new HashSet<>();
        try(DelimFileReader reader = new DelimFileReader(tsvFile))
        {
            for(DelimFileReader.Row row : reader)
            {
                String contig = row.getString(Columns.ref_contig);
                if(!contigs.add(contig))
                {
                    throw new UserInputError("Duplicate contig: " + contig);
                }
                OncologyGroup oncologyGroup = new OncologyGroup(row.getString(Columns.oncology_group));
                VirusInfo virusInfo = new VirusInfo(contig, row.getString(Columns.virus_name), oncologyGroup);
                result.add(virusInfo);
            }
        }
        return result;
    }

    private enum Columns
    {
        ref_contig,
        virus_name,
        oncology_group
    }
}
