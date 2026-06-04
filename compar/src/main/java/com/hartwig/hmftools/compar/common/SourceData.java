package com.hartwig.hmftools.compar.common;

import java.util.Map;

import com.google.common.collect.Maps;
import com.hartwig.hmftools.patientdb.dao.DatabaseAccess;

public class SourceData
{
    public final SourceType Type;
    public final DatabaseAccess Database;
    public final FileSources Files;

    public final Map<String,String> SampleIdMapping; // maps from the common tumor SampleId to the tumor Sample Id for this source
    public final Map<String,String> ReferenceSampleIdMapping; // as above for the reference Id

    public SourceData(final SourceType type, final DatabaseAccess database, final FileSources files)
    {
        Type = type;
        Database = database;
        Files = files;

        SampleIdMapping = Maps.newHashMap();
        ReferenceSampleIdMapping = Maps.newHashMap();
    }

    public String configName() { return Type.toString().toLowerCase(); }

    public String toString() { return Type.toString(); }
}
