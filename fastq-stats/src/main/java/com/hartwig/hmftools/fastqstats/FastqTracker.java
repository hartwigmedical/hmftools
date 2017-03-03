package com.hartwig.hmftools.fastqstats;

import java.util.SortedMap;

public class FastqTracker {
    private SortedMap<TrackerKey, Tracker> trackers;

    public FastqTracker(SortedMap<TrackerKey, Tracker> trackers) {
        this.trackers = trackers;
    }

    public void addValue(TrackerKey[] keys, int v){
        for(TrackerKey key: keys) {
            trackers.get(key).addValue(v);
        }
    }

    public long get(TrackerKey key){
        return trackers.get(key).getCount();
    }

    public double getPercentage(TrackerKey key, TrackerKey total){
        return trackers.get(key).getCount() * 100.0 / trackers.get(total).getCount();
    }

    public void putIfAbsent(TrackerKey key, Tracker t){
        trackers.putIfAbsent(key, t);
    }

    public String toString(){
        StringBuilder s = new StringBuilder();
        for(TrackerKey key: trackers.keySet()){
            s.append(key).append(": ").append(trackers.get(key).getCount()).append("\n");
        }
        return s.toString();
    }
}
