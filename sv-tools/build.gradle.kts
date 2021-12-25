
plugins {
    id("com.hartwig.java-conventions")
}

description = "HMF Tools - SV Tools"

dependencies {
    implementation(project(":hmf-common"))
    implementation(project(":patient-db"))
    annotationProcessor(libs.immutables.value)
    compileOnly(libs.immutables.value)
    testImplementation(libs.junit)
}

shadowJar {
    mainClass.set("com.hartwig.hmftools.svtools.simulation.ShatteringSim")
}
