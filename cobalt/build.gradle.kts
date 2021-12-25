
plugins {
    id("com.hartwig.java-conventions")
}

description = "HMF Tools - COBALT"

dependencies {
    implementation(project(":hmf-common"))
    testImplementation(libs.junit)
}

shadowJar {
    mainClass.set("com.hartwig.hmftools.cobalt.CobaltApplication")
}
