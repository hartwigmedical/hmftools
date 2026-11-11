package com.hartwig.hmftools.viridian.reference;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import com.hartwig.hmftools.common.utils.file.DelimFileReader;
import com.hartwig.hmftools.viridian.common.UserInputError;

// Per-virus-genome reference data.
public record VirusGenomeInfo(
        String contigName,
        Virus virus
)
{
    public static List<VirusGenomeInfo> load(String tsvFile)
    {
        List<VirusGenomeInfo> result = new ArrayList<>();
        Set<String> contigs = new HashSet<>();
        Map<String, Virus> virusesByName = new HashMap<>();
        try(DelimFileReader reader = new DelimFileReader(tsvFile))
        {
            for(DelimFileReader.Row row : reader)
            {
                String contig = row.getString(Column.contig_name);
                if(!contigs.add(contig))
                {
                    throw new UserInputError("Duplicate contig_name: " + contig);
                }
                VirusGroup virusGroup = new VirusGroup(row.getString(Column.virus_group));
                Virus virus = new Virus(row.getString(Column.virus_name), virusGroup);
                Virus existing = virusesByName.putIfAbsent(virus.name(), virus);
                if(existing != null && !existing.equals(virus))
                {
                    throw new UserInputError("Virus in more than one virus group: " + virus.name());
                }
                result.add(new VirusGenomeInfo(contig, virus));
            }
        }
        return result;
    }

    private enum Column
    {
        contig_name,
        virus_name,
        virus_group
    }
}
