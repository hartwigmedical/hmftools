package com.hartwig.hmftools.virusinterpreter.taxonomy;

import static com.hartwig.hmftools.virusinterpreter.VirusInterpreterApplication.VI_LOGGER;

import java.util.Map;

import com.google.common.annotations.VisibleForTesting;

import org.apache.logging.log4j.util.Strings;

public class TaxonomyDb
{
    private final Map<Integer, String> taxidToNameMap;

    public TaxonomyDb(final Map<Integer, String> taxidToNameMap)
    {
        this.taxidToNameMap = taxidToNameMap;
    }

    public String lookupName(int taxid)
    {
        boolean isMappedTaxid = taxidExists(taxid);

        if(!isMappedTaxid)
        {
            VI_LOGGER.warn("Could not match taxid '{}' to a name in taxonomy DB", taxid);
        }

        return isMappedTaxid ? taxidToNameMap.get(taxid) : Strings.EMPTY;
    }

    @VisibleForTesting
    boolean taxidExists(int taxid)
    {
        return taxidToNameMap.containsKey(taxid);
    }

    @VisibleForTesting
    int count()
    {
        return taxidToNameMap.keySet().size();
    }

}
