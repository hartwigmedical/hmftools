package com.hartwig.hmftools.common.amber;

import java.util.Collection;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Maps;

public class AmberPatientFactory
{
    private int mNextPatientId;
    private final Collection<AmberMapping> mMappings;
    private final Map<String, Integer> mPatientMap;

    public AmberPatientFactory(final List<AmberPatient> currentPatients, final Collection<AmberMapping> mappings)
    {
        mPatientMap = Maps.newHashMap();

        mMappings = mappings;
        for(AmberPatient patient : currentPatients)
        {
            mPatientMap.put(patient.sample(), patient.patientId());
        }

        mNextPatientId = currentPatients.stream().mapToInt(AmberPatient::patientId).max().orElse(0) + 1;
    }

    public AmberPatient createPatient(AmberSample sample)
    {
        final String sampleId = sample.sampleId();
        if(mPatientMap.containsKey(sampleId))
        {
            return create(sampleId, mPatientMap.get(sampleId));
        }

        for(AmberMapping mapping : mMappings)
        {
            if(mapping.firstSample().equals(sampleId) && mPatientMap.containsKey(mapping.secondSample()))
            {
                return create(sampleId, mPatientMap.get(mapping.secondSample()));
            }

            if(mapping.secondSample().equals(sampleId) && mPatientMap.containsKey(mapping.firstSample()))
            {
                return create(sampleId, mPatientMap.get(mapping.firstSample()));
            }
        }

        mPatientMap.put(sampleId, mNextPatientId);
        AmberPatient result = create(sampleId, mNextPatientId);
        mNextPatientId++;
        return result;
    }

    private static AmberPatient create(String sampleId, int patientId)
    {
        return ImmutableAmberPatient.builder()
                .patientId(patientId)
                .sample(sampleId)
                .build();
    }

}
