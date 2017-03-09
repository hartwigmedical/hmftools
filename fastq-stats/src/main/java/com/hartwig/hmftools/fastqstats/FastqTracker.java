package com.hartwig.hmftools.fastqstats;

import java.util.Map;
import java.util.TreeMap;

public class FastqTracker {
    private FastqData flowcell;
    private TreeMap<String, FastqData> lanes;
    private TreeMap<String, FastqData> samples;
    private FastqData undetermined;

    public FastqTracker(){
        flowcell = new FastqData(0, 0);
        undetermined = new FastqData(0, 0);
        lanes = new TreeMap<>();
        samples = new TreeMap<>();
    }

    public void addToFlowcell(FastqData data){
        flowcell = flowcell.add(data);
    }

    public void addToUndetermined(FastqData data) {undetermined = undetermined.add(data);}

    public void addToLane(String lane, FastqData data){
        if(!lanes.containsKey(lane)){
            lanes.put(lane, data);
        }
        else {
            FastqData current = lanes.get(lane);
            lanes.put(lane, current.add(data));
        }
    }

    public void addToSample(String sample, FastqData data){
        if(!samples.containsKey(sample)){
            samples.put(sample, data);
        }
        else {
            FastqData current = samples.get(sample);
            samples.put(sample, current.add(data));
        }
    }

    public FastqData getFlowcellData(){
        return flowcell;
    }

    public FastqData getUndeterminedData(){
        return undetermined;
    }

    public FastqData getLaneData(String lane){
        return lanes.get(lane);
    }

    public FastqData getSampleData(String sample){
        return samples.get(sample);
    }

    public Map<String, FastqData> getLanes(){
        return lanes;
    }

    public Map<String, FastqData> getSamples(){
        return samples;
    }
}
