package com.hartwig.hmftools.common.amber;

import java.util.Collection;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Maps;

public class AmberPatientFactory
{
    private int nextPatientId;
    private final Collection<AmberMapping> mappings;
    private final Map<String, Integer> patientMap = Maps.newHashMap();

    public AmberPatientFactory(final List<AmberPatient> currentPatients, final Collection<AmberMapping> mappings)
    {
        this.mappings = mappings;
        for(AmberPatient patient : currentPatients)
        {
            patientMap.put(patient.sample(), patient.patientId());
        }

        nextPatientId = currentPatients.stream().mapToInt(AmberPatient::patientId).max().orElse(0) + 1;
    }

    public AmberPatient createPatient(AmberSample sample)
    {
        final String sampleId = sample.sampleId();
        if(patientMap.containsKey(sampleId))
        {
            return create(sampleId, patientMap.get(sampleId));
        }

        for(AmberMapping mapping : mappings)
        {
            if(mapping.firstSample().equals(sampleId) && patientMap.containsKey(mapping.secondSample()))
            {
                return create(sampleId, patientMap.get(mapping.secondSample()));
            }

            if(mapping.secondSample().equals(sampleId) && patientMap.containsKey(mapping.firstSample()))
            {
                return create(sampleId, patientMap.get(mapping.firstSample()));
            }
        }

        patientMap.put(sampleId, nextPatientId);
        AmberPatient result = create(sampleId, nextPatientId);
        nextPatientId++;
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
