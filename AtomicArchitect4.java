

/*
Copyright (c) 2026 Sam78887
Licensed under SNSL v1.0
Created by: Sam78887
Non-Commercial Use Only
*/



/**

Atomic Lab Architect JAVA - ALAJ
 * Atomic Architect
 * A Java-based chemical analysis tool and molecular dynamics simulator.
 * * Features
 * * Electron Configuration
 * Automatically generates subshell notation (e.g., 1s^2 2s^2) for elements 1 through 118.
 * * VSEPR Theory
 * Predicts molecular geometry, bond angles, and steric numbers to determine the 3D shape.
 * * MD Simulation
 * Runs a digital lab simulation using a Lennard-Jones potential to model particle interactions.
 * * Physical State Prediction
 * Estimates boiling points and predicts whether a compound will be a solid, liquid, or gas.
 */
import java.util.*;

public class AtomicArchitect4  {

    static class Atom {
        String name, symbol, config;
        int valenceElectrons, valency, atomicNum;
        double weight, en, radius;

        static final String[] NAMES = {"", "Hydrogen", "Helium", "Lithium", "Beryllium", "Boron", "Carbon", "Nitrogen", "Oxygen", "Fluorine", "Neon", "Sodium", "Magnesium", "Aluminum", "Silicon", "Phosphorus", "Sulfur", "Chlorine", "Argon", "Potassium", "Calcium", "Scandium", "Titanium", "Vanadium", "Chromium", "Manganese", "Iron", "Cobalt", "Nickel", "Copper", "Zinc", "Gallium", "Germanium", "Assenic", "Selenium", "Bromine", "Krypton", "Rubidium", "Strontium", "Yttrium", "Zirconium", "Niobium", "Molybdenum", "Technetium", "Ruthenium", "Rhodium", "Palladium", "Silver", "Cadmium", "Indium", "Tin", "Antimony", "Tellurium", "Iodine", "Xenon", "Cesium", "Barium", "Lanthanum", "Cerium", "Praseodymium", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th", "Pa", "U", "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm", "Md", "No", "Lr", "Rf", "Db", "Sg", "Bh", "Hs", "Mt", "Ds", "Rg", "Cn", "Nh", "Fl", "Mc", "Lv", "Ts", "Og"};
        static final String[] SYMBOLS = {"", "H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne", "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar", "K", "Ca", "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr", "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn", "Sb", "Te", "I", "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th", "Pa", "U", "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm", "Md", "No", "Lr", "Rf", "Db", "Sg", "Bh", "Hs", "Mt", "Ds", "Rg", "Cn", "Nh", "Fl", "Mc", "Lv", "Ts", "Og"};
        static final double[] WEIGHTS = new double[119];
        static final double[] EN_SCALE = new double[119];
        static final double[] RADII = new double[119];

static {
    // 1. Initialize with defaults to prevent nulls Created by: Sam78887
    Arrays.fill(WEIGHTS, 0.0);
    Arrays.fill(EN_SCALE, 0.0);
    Arrays.fill(RADII, 0.0);

    // 2. Full Atomic Weights (Standard Atomic Masses) Created by: Sam78887
    double[] allWeights = {
        0, 1.008, 4.0026, 6.94, 9.0122, 10.81, 12.011, 14.007, 15.999, 18.998, 20.18, // 0-10
        22.99, 24.305, 26.982, 28.085, 30.974, 32.06, 35.45, 39.948, 39.098, 40.078, // 11-20
        44.956, 47.867, 50.942, 51.996, 54.938, 55.845, 58.933, 58.693, 63.546, 65.38, // 21-30
        69.723, 72.63, 74.922, 78.971, 79.904, 83.798, 85.468, 87.62, 88.906, 91.224, // 31-40
        92.906, 95.95, 98, 101.07, 102.91, 106.42, 107.87, 112.41, 114.82, 118.71, // 41-50
        121.76, 127.6, 126.9, 131.29, 132.91, 137.33, 138.91, 140.12, 140.91, 144.24, // 51-60
        145, 150.36, 151.96, 157.25, 158.93, 162.5, 164.93, 167.26, 168.93, 173.05, // 61-70
        174.97, 178.49, 180.95, 183.84, 186.21, 190.23, 192.22, 195.08, 196.97, 200.59, // 71-80
        204.38, 207.2, 208.98, 209, 210, 222, 223, 226, 227, 232.04, 231.04, 238.03, // 81-92
        237, 244, 243, 247, 247, 251, 252, 257, 258, 259, 262, 267, 268, 271, 270, // 93-108
        270, 278, 281, 281, 285, 286, 289, 289, 293, 294, 294 // 109-118
    };
    System.arraycopy(allWeights, 0, WEIGHTS, 0, allWeights.length);

    // 3. Pauling Electronegativity (EN_SCALE)
    // Common elements manually set, others follow a calculated periodic trend Created by: Sam78887
    double[] en = {
        0, 2.20, 0, 0.98, 1.57, 2.04, 2.55, 3.04, 3.44, 3.98, 0, // 0-10
        0.93, 1.31, 1.61, 1.90, 2.19, 2.58, 3.16, 0, 0.82, 1.00 // 11-20
    };
    System.arraycopy(en, 0, EN_SCALE, 0, en.length);
    // Fill remaining EN based on general trends (Group 17 high, Group 1 low)
    for(int i=21; i<=118; i++) {
        if (EN_SCALE[i] == 0) {
            if (i >= 57 && i <= 71) EN_SCALE[i] = 1.1; // Lanthanides
            else if (i >= 89 && i <= 103) EN_SCALE[i] = 1.3; // Actinides
            else EN_SCALE[i] = 1.5 + (0.01 * (i % 18)); // Rough heuristic for heavy metals
        }
    }
    EN_SCALE[47]=1.93; EN_SCALE[79]=2.54; EN_SCALE[17]=3.16; EN_SCALE[35]=2.96; EN_SCALE[53]=2.66;

    // 4. Atomic Radii (pm)
    // Simplified Radii trends: Decreases across period, increases down group
    for (int i = 1; i <= 118; i++) {
        int shell = (i <= 2) ? 1 : (i <= 10) ? 2 : (i <= 18) ? 3 : (i <= 36) ? 4 : (i <= 54) ? 5 : (i <= 86) ? 6 : 7;
        RADII[i] = 30 + (shell * 25) - ((i % 10) * 2);
    }
    // Specific Overrides for accuracy Created by: Sam78887
    RADII[1]=37; RADII[6]=77; RADII[7]=75; RADII[8]=73; RADII[9]=71; 
    RADII[17]=99; RADII[47]=144; RADII[79]=144; RADII[80]=151;
}

        Atom(int num) {
            this.atomicNum = (num >= 1 && num <= 118) ? num : 0;
            this.name = NAMES[this.atomicNum];
            this.symbol = SYMBOLS[this.atomicNum];
            this.weight = WEIGHTS[this.atomicNum];
            this.en = EN_SCALE[this.atomicNum];
            this.radius = RADII[this.atomicNum];
            generateConfigAndValence(this.atomicNum);
        }

private void generateConfigAndValence(int n) {
    String[] subshells = {"1s", "2s", "2p", "3s", "3p", "4s", "3d", "4p", "5s", "4d", "5p", "6s", "4f", "5d", "6p", "7s", "5f", "6d", "7p"};
    int[] caps = {2, 2, 6, 2, 6, 2, 10, 6, 2, 10, 6, 2, 14, 10, 6, 2, 14, 10, 6};
    StringBuilder sb = new StringBuilder();
    int temp = n;

    // 1. Standard Aufbau Filling Loop Created by: Sam78887
    for (int i = 0; i < subshells.length && temp > 0; i++) {
        int fill = Math.min(temp, caps[i]);
        sb.append(subshells[i]).append("^").append(fill).append(" ");
        temp -= fill;
    }
    this.config = sb.toString().trim();

    // 2. Exception Handler (Stability of Half/Full d-shells)
    // These override the string generated by the loop above Created by: Sam78887
    if (n == 24) this.config = "1s^2 2s^2 2p^6 3s^2 3p^6 4s^1 3d^5";      // Chromium (Cr)
    else if (n == 29) this.config = "1s^2 2s^2 2p^6 3s^2 3p^6 4s^1 3d^10";     // Copper (Cu)
    else if (n == 42) this.config = "1s^2 2s^2 2p^6 3s^2 3p^6 4s^2 3d^10 4p^6 5s^1 4d^5";  // Molybdenum (Mo)
    else if (n == 47) this.config = "1s^2 2s^2 2p^6 3s^2 3p^6 4s^2 3d^10 4p^6 5s^1 4d^10"; // Silver (Ag)
    else if (n == 79) this.config = "[Xe] 4f^14 5d^10 6s^1"; // Gold (Au)

    // 3. Valence and Valency Logic
    if (n == 2 || n == 10 || n == 18 || n == 36 || n == 54 || n == 86 || n == 118) {
        this.valency = 0; 
        this.valenceElectrons = (n == 2) ? 2 : 8; // Correction: Helium has 2 valence electrons
    } else {
        int mod = n % 18;
        this.valenceElectrons = (mod == 0) ? 8 : (mod > 10 ? mod - 10 : (mod <= 2 ? mod : (mod <= 12 ? 2 : mod - 10)));
        
        // Manual Overrides for specific group behaviors Created by: Sam78887
        if (n == 1) valenceElectrons = 1; 
        if (n == 6) valenceElectrons = 4;
        if (n == 7) valenceElectrons = 5;
        if (n == 8) valenceElectrons = 6;
        if (n == 9) valenceElectrons = 7;
        
        // Correcting Valence Electrons for the transition metals with exceptions Created by: Sam78887
        if (n == 24 || n == 29 || n == 47 || n == 79) valenceElectrons = 1;

        this.valency = (valenceElectrons <= 4) ? valenceElectrons : 8 - valenceElectrons;
    }
}

    }
    static class Particle {
        double x, y, z, vx, vy, vz, fx, fy, fz;
        Particle(double x, double y, double z) { this.x = x; this.y = y; this.z = z; }
    }
    public static void main(String[] args) {
        Scanner sc = new Scanner(System.in);
        String choice = "y";
        while (choice.equalsIgnoreCase("y")) {
            // 1. Initialization
            Map<String, Integer> counts = new LinkedHashMap<>();
            double totalWeight = 0;
            int totalValencySum = 0;
            int maxVal = -1;
            Atom central = null, secondary = null;
            int totalAtoms = 0;

            System.out.print("\nHow many unique element types? ");
            int typeCount = sc.nextInt();

            // 2. Element Input Loop
            for (int i = 0; i < typeCount; i++) {
                System.out.print("Enter Atomic Number (1-118): ");
                int z = sc.nextInt();
                System.out.print("Enter Quantity: ");
                int qty = sc.nextInt();
                
                Atom cur = new Atom(z);
                System.out.println("-> " + cur.name + " [" + cur.config + "]");
                
                totalWeight += (cur.weight * qty);
                totalValencySum += (cur.valency * qty);
                totalAtoms += qty;
                counts.put(cur.symbol, qty);

                // Central Atom Logic Created by: Sam78887
                if (cur.valency > maxVal || (cur.valency == maxVal && qty == 1)) {
                    maxVal = cur.valency;
                    secondary = central;
                    central = cur;
                } else if (secondary == null && !cur.symbol.equals(central.symbol)) {
                    secondary = cur;
                }
            }

            // 3. Phase & Discovery Analysis (Always Runs)
            predictPhysicalState(counts, totalWeight);

            System.out.println("\n--- Discovery Scanner Analysis ---");
            if (totalAtoms == 1) {
                System.out.println("Status: [PURE ELEMENTAL LATTICE / MONATOMIC]");
            } else {
                boolean valCheck = (totalValencySum % 2 == 0);
                boolean saturated = true;
                if (counts.containsKey("C")) {
                    int others = totalAtoms - counts.get("C");
                    if (others < (counts.get("C") / 2 + 1) && !counts.containsKey("O")) saturated = false;
                }
                
                if (!valCheck || totalValencySum < (2 * maxVal)) {
                    System.out.println("Status: [IMPOSSIBLE / UNSTABLE RADICAL]");
                } else if (!saturated) {
                    System.out.println("Status: [HIGHLY REACTIVE / THEORETICAL]");
                } else {
                    System.out.println("Status: [STABLE COMPOUND]");
                }
            }

            // 4. Composition & Structural Details Created by: Sam78887
            System.out.println("\n--- Mass Composition ---");
            for (String s : counts.keySet()) {
                double atomW = 0;
                for (int i=1; i<119; i++) if (Atom.SYMBOLS[i].equals(s)) atomW = Atom.WEIGHTS[i];
                System.out.printf("%s: %.2f%%\n", s, (atomW * counts.get(s) / totalWeight) * 100);
            }

            System.out.println("\n--- Structural Details ---");
            System.out.println("Molecular Formula: " + generateFormula(counts));
            System.out.println("SMILES Notation: " + generateSMILES(central, counts));
            System.out.printf("Total Molecular Weight: %.5f g/mol\n", totalWeight);

            // 5. Reaction Test (Independent)
            System.out.print("\nAdd a second molecule for a reaction test? (y/n): ");
            if (sc.next().equalsIgnoreCase("y")) {
                Map<String, Integer> countsB = new HashMap<>();
                System.out.print("Second Molecule - Unique element types: ");
                int typesB = sc.nextInt();
                for (int j = 0; j < typesB; j++) {
                    System.out.print("Atomic Num (1-118): ");
                    int zB = sc.nextInt();
                    System.out.print("Quantity: ");
                    int qB = sc.nextInt();
                    countsB.put(Atom.SYMBOLS[zB], qB);
                }
                predictReaction(counts, countsB);
            }

            // 6. Bond Analysis Block (Only for Compounds) Created by: Sam78887
            if (secondary != null) {
                System.out.println("\n--- Bond Analysis ---");
                double diff = Math.abs(central.en - secondary.en);
                String bType = (diff >= 1.7) ? "Ionic" : (diff >= 0.5) ? "Polar Covalent" : "Non-Polar Covalent";
                System.out.println("Primary Bond Type: " + bType + " (Delta EN: " + String.format("%.2f", diff) + ")");
                
                double bondLength = (central.radius + secondary.radius) - 9 * diff;
                System.out.printf("Predicted Bond Length (%s-%s): %.2f pm\n", central.symbol, secondary.symbol, bondLength);
                
                // Spectroscopy logic
                double muAmu = (central.weight * secondary.weight) / (central.weight + secondary.weight);
                double muKg = muAmu * 1.660539e-27;
                double k = 500.0 * (1.0 + Math.abs(central.en - secondary.en) / 2.0);
                double freqHz = (1.0 / (2.0 * Math.PI)) * Math.sqrt(k / muKg);
                System.out.printf("Vibration Frequency: %.2f THz\n", freqHz / 1.0e12);

                if (counts.containsKey("O") && counts.containsKey("H")) {
                    System.out.println("Intermolecular Force: Strong Hydrogen Bonding detected.");
                }
            }

            // 7. Multi-Center Geometry (Accessible to all) Created by: Sam78887
            if (totalAtoms >= 1) {
                System.out.println("\n--- Multi-Center Geometry Analysis ---");
                if (totalAtoms == 1) {
                    System.out.println("Geometry: Point / Monatomic Particle");
                } else {
                    if (counts.getOrDefault("C", 0) > 0) System.out.println("At Carbon Center: Tetrahedral (109.5°)");
                    if (counts.containsKey("O") && counts.containsKey("H")) {
                        System.out.print("At Oxygen (-OH) Center: ");
                        predictGeometry(4, 2); // Steric 4, LP 2 = Bent Created by: Sam78887
                    }
                }
            }

            // 8. Digital Lab Simulation (Accessible to all)
            System.out.print("\nRun Digital Lab Simulation (MD) for this molecule? (y/n): ");
            if (sc.next().equalsIgnoreCase("y")) {
                runLabSimulation(counts, totalWeight);
            }
            
            System.out.print("\nAnalyze another? (y/n): ");
            choice = sc.next();
        }

        }


    public static void runLabSimulation(Map<String, Integer> counts, double molWeight) {
        List<Particle> particles = new ArrayList<>();
        double sigma = 3.5, epsilon = 0.2, targetTemp = 200.0;
        String molName = "Generic Compound";

        // Inside runLabSimulation Created by: Sam78887
int cCount = counts.getOrDefault("C", 0);
boolean hBond = counts.containsKey("H") && (counts.containsKey("O") || counts.containsKey("N"));

// Sync the Lab Thermostat with the Physical State Prediction
double bpCelsius = -170.0 + (cCount * 32.0) + (molWeight * 0.2) + (hBond ? 160.0 : 0.0);
if (cCount == 0 && counts.getOrDefault("O", 0) == 1) bpCelsius = 100.0; // Water override

targetTemp = bpCelsius + 273.15; // Convert to Kelvin for the MD loop

// Calibration: Adjust "Stickiness" (Epsilon) based on chain length
// Created by: Sam78887 Created by: Sam78887 Created by: Sam78887
// This prevents the "Freezing" effect in the small 40-molecule box
epsilon = hBond ? 0.55 : (cCount > 4 ? 0.22 : 0.15);


        // Heuristic Intelligence Layer
        // Created by: Sam78887
        boolean hBonding = counts.containsKey("H") && (counts.containsKey("N") || counts.containsKey("O") || counts.containsKey("F"));
        if (counts.getOrDefault("C",0)==1 && counts.getOrDefault("H",0)==4) { epsilon=0.148; sigma=3.73; targetTemp=111.6; molName="Methane (CH4)"; }
        else if (counts.getOrDefault("O",0)==1 && counts.getOrDefault("H",0)==2) { epsilon=0.650; sigma=3.16; targetTemp=373.15; molName="Water (H2O)"; }
        else {

            
            // 1. Get structural counts for the formula
            
            int hCount = counts.getOrDefault("H", 0);
            int oCount = counts.getOrDefault("O", 0);
            int nCount = counts.getOrDefault("N", 0);
            int sCount = counts.getOrDefault("S", 0);

            // 2. Linear Surface Area Logic: +22.0K per Carbon fixes the "Heptanol Gap"
            // Base value (40.0) + Carbon chain influence + heavy atom mass influence
            // Created by: Sam78887
            double surfaceAreaBase = (cCount * 22.0) + (molWeight * 0.4) + 40.0;

            // 3. Polarity Logic: Differentiates OH (Alcohol) vs NH (Amine) vs SH (Thiol)
            double polarityBonus = 0;
            if (hCount > 0) {
                if (oCount > 0) polarityBonus = 165.0;      // Strong H-Bonding
                else if (nCount > 0) polarityBonus = 85.0;  // Medium H-Bonding
                else if (sCount > 0) polarityBonus = 20.0;  // Very Weak (Sulfur)
            }

            // 4. Final Target Calculation
            targetTemp = surfaceAreaBase + polarityBonus + 200.0;
            
            // Adjust Epsilon (Stickiness) based on polarity
            // Created by: Sam78887
            epsilon = (oCount > 0 || nCount > 0) ? 0.55 : 0.25; 
            
            molName = generateFormula(counts);
        }


        System.out.println("\n--- STABILIZED DIGITAL LAB: " + molName + " ---");
        int numMols = 40; double dt = 0.5;
        int side = (int) Math.ceil(Math.pow(numMols, 1.0/3.0));
        for (int i = 0; i < numMols; i++) particles.add(new Particle((i%side)*6.0, ((i/side)%side)*6.0, (i/(side*side))*6.0));

        for (int step = 0; step < 1001; step++) {
            for (Particle p : particles) { p.fx = p.fy = p.fz = 0; }
            for (int i = 0; i < particles.size(); i++) {
                for (int j = i + 1; j < particles.size(); j++) {
                    Particle p1 = particles.get(i), p2 = particles.get(j);
                    double dx = p1.x-p2.x, dy = p1.y-p2.y, dz = p1.z-p2.z;
                    double r2 = dx*dx + dy*dy + dz*dz + 0.5;
                    double r = Math.sqrt(r2);
                    if (r < 12.0) {
                        double r6 = Math.pow(sigma / r, 6);
                        double fS = (48 * epsilon / r2) * (r6 * r6 - 0.5 * r6);
                        if (fS > 80) fS = 80;
                        p1.fx += fS * dx; p1.fy += fS * dy; p1.fz += fS * dz;
                        p2.fx -= fS * dx; p2.fy -= fS * dy; p2.fz -= fS * dz;
                    }
                }
            }
            double currentKE = 0;
            for (Particle p : particles) {
                double acc = 0.0004184 / molWeight;
                p.vx += p.fx * acc * dt; p.vy += p.fy * acc * dt; p.vz += p.fz * acc * dt;
                p.vx = Math.max(-1.5, Math.min(1.5, p.vx));
                p.x += p.vx * dt; p.y += p.vy * dt; p.z += p.vz * dt;
                currentKE += (p.vx*p.vx + p.vy*p.vy + p.vz*p.vz);
            }
            double currentT = (currentKE * molWeight) / (3 * numMols * 0.008314);
            if (step % 10 == 0) {
                double scale = Math.sqrt(1.0 + 0.1 * (targetTemp / (currentT + 0.1) - 1.0));
                for (Particle p : particles) { p.vx *= scale; p.vy *= scale; p.vz *= scale; }
            }
            if (step % 200 == 0) System.out.printf("Step %d | Temp: %.2f K | Status: %s\n", step, currentT, (currentT >= targetTemp - 5 ? "STABILIZED" : "HEATING"));
        }
    }
public static void predictPhysicalState(Map<String, Integer> counts, double weight) {
    // 1. Setup Totals
    int totalAtoms = 0;
    for (int qty : counts.values()) totalAtoms += qty;
    
    int cCount = counts.getOrDefault("C", 0);
    int hCount = counts.getOrDefault("H", 0);
    int oCount = counts.getOrDefault("O", 0);
    int nCount = counts.getOrDefault("N", 0);
    String mainSymbol = (counts.size() == 1) ? counts.keySet().iterator().next() : "";

    // 2. Identity Flags
    boolean isNobleGas = "He Ne Ar Kr Xe Rn".contains(mainSymbol) && counts.size() == 1;
    boolean isMetal = "Li Na K Mg Al Fe Cu Zn Ag Au Pb Sn W".contains(mainSymbol) && counts.size() == 1;

    double bpCelsius;
    String predictedPhase;

    // 3. Logic Engine (Hierarchy of Exceptions)
    if (mainSymbol.equals("Hg") && totalAtoms == 1) {
        bpCelsius = 356.7; // Mercury specific
        predictedPhase = "Liquid (at 25°C)";
    } 
    else if (isMetal) {
        // Metallic bonding is extremely strong (Sea of Electrons)
        bpCelsius = 1500.0 + (weight * 5.0); 
        if (mainSymbol.equals("W")) bpCelsius = 5555.0; // Tungsten Exception
        predictedPhase = "Solid (Metallic Lattice)";
    } 
    else if (isNobleGas) {
        bpCelsius = -272.0 + (weight * 0.4);
        predictedPhase = "Gas (at 25°C)";
    } 
    else if (cCount == 0 && oCount == 1 && hCount == 2 && totalAtoms == 3) {
        bpCelsius = 100.0; // Water override
        predictedPhase = "Liquid (at 25°C)";
    } 
    //noble gas condition
    // Created by: Sam78887

    else if (isNobleGas) {
        // Refined Heuristic for London Dispersion Forces
        // Works for He (-269C) up to Rn (-61C)
        bpCelsius = -272.0 + (weight * 1.26);
        
        // Specific Overrides for precision
        if (mainSymbol.equals("He")) bpCelsius = -268.9;
        if (mainSymbol.equals("Ne")) bpCelsius = -246.1;
        
        predictedPhase = (bpCelsius > 25) ? "Liquid" : "Gas (at 25°C)";
    }
    else {
        // Organic / Covalent Logic
        boolean hBond = hCount > 0 && (oCount > 0 || nCount > 0);
        double baseBP = -170.0; 
        double chainEffect = (cCount * 32.0); 
        double massEffect = (weight * 0.22);
        double polarityEffect = hBond ? 160.0 : 0.0;
        
        bpCelsius = baseBP + chainEffect + massEffect + polarityEffect;

        // Phase Check for Covalent
        if (bpCelsius < 25) predictedPhase = "Gas (at 25°C)";
        else if (cCount > 12 || weight > 300) predictedPhase = "Solid (Organic/Wax)";
        else predictedPhase = "Liquid (at 25°C)";
    }

    // 4. Output
    System.out.println("\n--- Physical State Prediction ---");
    System.out.printf("Estimated Boiling Point: ~%.1f°C\n", bpCelsius);
    System.out.println("Predicted Phase: " + predictedPhase);
}
public static String generateFormula(Map<String, Integer> counts) {
    // Specialized format for Alcohols (e.g., C2H5OH instead of C2H6O)
    if (counts.getOrDefault("C", 0) > 0 && counts.getOrDefault("O", 0) > 0 && counts.getOrDefault("H", 0) > 1) {
        int c = counts.get("C");
        int h = counts.get("H");
        return String.format("C%sH%sOH", (c > 1 ? c : ""), (h - 1 > 1 ? (h - 1) : ""));
    }

    // Hill System: Carbon first, then Hydrogen, then Alphabetical
    StringBuilder sb = new StringBuilder();
    List<String> sorted = new ArrayList<>(counts.keySet());
    Collections.sort(sorted, (a, b) -> {
        if (a.equals("C")) return -1;
        if (b.equals("C")) return 1;
        if (a.equals("H")) return -1;
        if (b.equals("H")) return 1;
        return a.compareTo(b);
    });

    for (String s : sorted) {
        int qty = counts.get(s);
        sb.append(s).append(qty > 1 ? qty : "");
    }
    return sb.toString();
}



public static String generateSMILES(Atom central, Map<String, Integer> counts) {
    StringBuilder sm = new StringBuilder();

    // 1. Handle Organic Chains (C-C-C...)
    if (counts.containsKey("C")) {
        for (int i = 0; i < counts.get("C"); i++) sm.append("C");
    } else {
        // For Inorganics, start with the Central Atom
        sm.append(central.symbol);
    }

    // 2. Add Functional Groups (O, N, S, P, etc.)
    // We skip Hydrogen in SMILES because it's usually "implicit"
    // Created by: Sam78887
    for (Map.Entry<String, Integer> entry : counts.entrySet()) {
        String sym = entry.getKey();
        if (!sym.equals("C") && !sym.equals("H") && !sym.equals(central.symbol)) {
            sm.append("(").append(sym).append(")");
        }
    }
    
    // 3. Special case for Alcohols: If Oxygen is present, ensure it shows as the -OH group
    if (counts.containsKey("O") && !sm.toString().contains("O")) {
        sm.append("O");
    }

    return sm.toString();
}


    // Updated predictGeometry method:
public static void predictGeometry(int sn, int lp) {
    if (sn == 2) System.out.println("Linear (180°)");
    else if (sn == 3) System.out.println(lp == 0 ? "Trigonal Planar (120°)" : "Bent (~118°)");
    else if (sn == 4) {
        if (lp == 0) System.out.println("Tetrahedral (109.5°)");
        else if (lp == 1) System.out.println("Trigonal Pyramidal (~107°)");
        else if (lp == 2) System.out.println("Bent (~104.5°) due to Lone Pair Repulsion");
    } else {
        System.out.println("Complex/Advanced Geometry");
    }
}


public static void predictReaction(Map<String, Integer> molA, Map<String, Integer> molB) {
    System.out.println("\n--- ⚗️ Reaction Prediction Lab ---");
    
// --- Molecule A Profile ---
    int cA = molA.getOrDefault("C", 0), hA = molA.getOrDefault("H", 0), oA = molA.getOrDefault("O", 0);
    int nA = molA.getOrDefault("N", 0), sA = molA.getOrDefault("S", 0), fA = molA.getOrDefault("F", 0);
    int clA = molA.getOrDefault("Cl", 0), pA = molA.getOrDefault("P", 0), iA = molA.getOrDefault("I", 0);
    int mgA = molA.getOrDefault("Mg", 0), alA = molA.getOrDefault("Al", 0), siA = molA.getOrDefault("Si", 0);

    // --- Molecule B Profile ---
    int cB = molB.getOrDefault("C", 0), hB = molB.getOrDefault("H", 0), oB = molB.getOrDefault("O", 0);
    int nB = molB.getOrDefault("N", 0), sB = molB.getOrDefault("S", 0), fB = molB.getOrDefault("F", 0);
    int clB = molB.getOrDefault("Cl", 0), pB = molB.getOrDefault("P", 0), iB = molB.getOrDefault("I", 0);
    int mgB = molB.getOrDefault("Mg", 0), alB = molB.getOrDefault("Al", 0), siB = molB.getOrDefault("Si", 0);

    // 1. Combustion Check (Fuel + Oxygen)
    if (cA > 0 && oB >= 2 && cB == 0) {
        System.out.println("Reaction: [COMBUSTION / OXIDATION]");
        double enthalpy = (cA * 393.5) + (hA * 141.8); 
        System.out.printf("Products: %dCO2 + %dH2O\n", cA, hA/2);
        System.out.printf("Energy Released: -%.2f kJ/mol (Highly Exothermic)\n", enthalpy);
    } 
    // Created by: Sam78887
    

    // 2. Fischer Esterification (Alcohol + Carboxylic Acid)
    else if(oA >= 1 && hA >= 1 && cB >= 1 && oB >= 2 && hB >= 1) {
        System.out.println("Reaction: [FISCHER ESTERIFICATION]");
        System.out.println("Mechanism: Nucleophilic Acyl Substitution");
        int totalC = cA + cB;
        int totalH = hA + hB - 2;
        System.out.printf("Product: Ester (C%dH%dO2) + Water (H2O)\n", totalC, totalH);
        
        // Aroma Database Logic
        // Created by: Sam78887
        if (totalC == 7) System.out.println("Aroma Detected: Banana (Isoamyl Acetate)");
        else if (totalC == 9) System.out.println("Aroma Detected: Orange (Octyl Acetate)");
        else System.out.println("Note: Likely produces a fruity or floral aroma.");
    }

    // 3. Acid-Base Neutralization (Amine + Acid)
    else if (nA >= 1 && clB >= 1) {
        System.out.println("Reaction: [ACID-BASE NEUTRALIZATION]");
        System.out.println("Mechanism: Proton Transfer (Brønsted–Lowry)");
        System.out.println("Product: Ammonium Chloride (NH4Cl) - White Salt/Smoke.");
    }

    // 4. Halogenation (Organic + Fluorine/Chlorine)
    else if (cA > 0 && (fB >= 2 || clB >= 2)) {
        System.out.println("Reaction: [HALOGENATION]");
        System.out.println("Mechanism: Free Radical Substitution");
        System.out.println("Status: Carbon-Halogen bonds forming. Potential Teflon precursor.");
    }

    // 5. Synthesis of Sulfides (Metal/H2 + Sulfur)
    else if (sB > 0 && (hA > 0 || molA.containsKey("Fe") || molA.containsKey("Cu"))) {
        System.out.println("Reaction: [SULFIDE FORMATION]");
        if (hA > 0) System.out.println("Product: Hydrogen Sulfide (H2S) - Warning: Toxic/Rotten Egg Odor.");
        else System.out.println("Product: Metallic Sulfide (Tarnish/Mineral Layer).");
    }

    // 6. Precipitation (Silver + Halogen) - Specific Test
    else if (molA.containsKey("Ag") && clB > 0) {
        System.out.println("Reaction: [PRECIPITATION]");
        System.out.println("Product: Silver Chloride (AgCl) - White insoluble solid.");
    }
    // 7. Alkali Metal + Water (The "Explosive" Test)
    else if ((molA.containsKey("Na") || molA.containsKey("K") || molA.containsKey("Li")) && hB == 2 && oB == 1) {
        System.out.println("Reaction: [VIOLENT ALKALI HYDRATION]");
        System.out.println("Mechanism: Single Displacement / Exothermic");
        System.out.println("Product: Metallic Hydroxide + Hydrogen Gas (H2).");
        System.out.println("Warning: Highly explosive; generates purple/orange flames.");
    }

    // 8. Stomach Antacid (Magnesium Hydroxide + HCl)
    else if (molA.containsKey("Mg") && oA >= 2 && clB >= 1) {
        System.out.println("Reaction: [ACID-BASE NEUTRALIZATION]");
        System.out.println("Product: Magnesium Chloride (Soluble Salt) + Water.");
        System.out.println("Medical Note: Relieves gastric hyperacidity.");
    }

    // 9. Thermite Reaction (Aluminum + Iron Oxide)
    else if (molA.containsKey("Al") && molB.containsKey("Fe") && oB >= 1) {
        System.out.println("Reaction: [THERMITE / REDOX]");
        System.out.println("Mechanism: Aluminum reducing Iron Oxide.");
        System.out.println("Product: Molten Iron + Aluminum Oxide.");
        System.out.println("Note: Used for welding railway tracks. Temp > 2500°C.");
    }

    // 10. Haber Process (Nitrogen + Hydrogen)
    else if (nA >= 2 && hB >= 2 && cA == 0 && cB == 0) {
        System.out.println("Reaction: [HABER-BOSCH SYNTHESIS]");
        System.out.println("Product: Ammonia (NH3).");
        System.out.println("Industrial Note: Essential for global fertilizer production.");
    }

    // 11. Rusting / Corrosion (Iron + Oxygen/Water)
    else if (molA.containsKey("Fe") && (oB >= 1 || (hB == 2 && oB == 1))) {
        System.out.println("Reaction: [OXIDATION / CORROSION]");
        System.out.println("Product: Iron(III) Oxide (Fe2O3) - Red Rust.");
        System.out.println("Structural Note: Causes significant material degradation over time.");
    }

    // 12. Carbonation (Water + CO2)
    else if (hA == 2 && oA == 1 && cB == 1 && oB == 2) {
        System.out.println("Reaction: [CARBONATION]");
        System.out.println("Product: Carbonic Acid (H2CO3).");
        System.out.println("Note: This is the 'fizz' in soft drinks and causes ocean acidification.");
    }

    // 13. Photosynthesis Simulation (CO2 + Water)
    else if (cA == 1 && oA == 2 && hB == 2 && oB == 1) {
        System.out.println("Reaction: [PHOTOSYNTHESIS (Light Required)]");
        System.out.println("Product: Glucose (C6H12O6) + Oxygen (O2).");
        System.out.println("Biosphere Note: Primary energy source for life on Earth.");
    }

    // 14. Copper Tarnish (Copper + Sulfur/Oxygen)
    else if (molA.containsKey("Cu") && (sB >= 1 || oB >= 1)) {
        System.out.println("Reaction: [OXIDATION / TARNISHING]");
        System.out.println("Product: Copper(II) Oxide or Copper Sulfide.");
        System.out.println("Observation: Color change to dull black or green (Patina).");
    }

    // 15. Etching (Silicon + Fluorine)
    else if (molA.containsKey("Si") && fB >= 2) {
        System.out.println("Reaction: [SEMICONDUCTOR ETCHING]");
        System.out.println("Product: Silicon Tetrafluoride (SiF4).");
        System.out.println("Tech Note: Critical process in microchip manufacturing.");
    }

    // 16. Bleaching (Sodium Hypochlorite + Stains)
    else if (molA.containsKey("Na") && clA == 1 && oA == 1 && cB >= 1) {
        System.out.println("Reaction: [OXIDATIVE BLEACHING]");
        System.out.println("Mechanism: Breaking of Chromophore double bonds.");
        System.out.println("Effect: Decolorization of organic pigments.");
    }

    // 17. Rocket Fuel (Hydrazine + Peroxide)
    else if (nA == 2 && hA == 4 && oB == 2 && hB == 2) {
        System.out.println("Reaction: [HYPERGOLIC COMBUSTION]");
        System.out.println("Product: Nitrogen Gas + Water Vapor.");
        System.out.println("Note: Spontaneous ignition used in spacecraft thrusters.");
    }

    // 18. Calcium Carbide (Carbide + Water)
    else if (molA.containsKey("Ca") && cA == 2 && hB == 2 && oB == 1) {
        System.out.println("Reaction: [ACETYLENE GENERATION]");
        System.out.println("Product: Acetylene Gas (C2H2) + Calcium Hydroxide.");
        System.out.println("Usage: Used in old miner lamps and industrial welding.");
    }
    // 19. Zinc + HCl (Standard Laboratory Hydrogen Generation)
    else if (molA.containsKey("Zn") && hB == 1 && clB == 1) {
        System.out.println("Reaction: [SINGLE DISPLACEMENT]");
        System.out.println("Product: Zinc Chloride (Salt) + Hydrogen Gas (H2 Bubbles).");
        System.out.println("Note: Standard method for 'Pop Test' hydrogen production.");
    }

    // 20. Lithium-Ion Discharge (Battery Simulation)
    else if (molA.containsKey("Li") && cA >= 6 && molB.containsKey("Co") && oB == 2) {
        System.out.println("Reaction: [ELECTROCHEMICAL DISCHARGE]");
        System.out.println("Mechanism: Lithium Intercalation.");
        System.out.println("Status: Battery powering device. Li moving to Cobalt Oxide cathode.");
    }

    // 21. Acid Rain Effect (Marble/Limestone + Acid)
    else if (molA.containsKey("Ca") && cA == 1 && oA == 3 && hB == 1 && clB == 1) {
        System.out.println("Reaction: [ACIDIC DEGRADATION]");
        System.out.println("Product: Calcium Chloride + Water + CO2.");
        System.out.println("Environmental Note: Simulates the erosion of marble statues by acid rain.");
    }

    // 22. Aqua Regia Component (Nitric + Hydrochloric)
    else if (nA == 1 && oA == 3 && hA == 1 && hB == 1 && clB == 1) {
        System.out.println("Reaction: [AQUA REGIA FORMATION]");
        System.out.println("Product: Nitrosyl Chloride (NOCl) + Free Chlorine.");
        System.out.println("Warning: Extremely corrosive mixture; can dissolve Gold and Platinum.");
    }

    // 23. Dehydration of Sugar (Sugar + Sulfuric Acid)
    else if (cA == 12 && hA == 22 && oA == 11 && sB == 1 && oB == 4) {
        System.out.println("Reaction: [DEHYDRATION / CARBON COLUMN]");
        System.out.println("Mechanism: Acid stripping water from carbohydrate.");
        System.out.println("Product: Pure Carbon (Black Foam) + Steam.");
        System.out.println("Observation: Rapid expansion of a black carbon 'snake'.");
    }

    // 24. Silver Tarnish Removal (Silver Sulfide + Aluminum)
    else if (molA.containsKey("Ag") && sA == 1 && molB.containsKey("Al")) {
        System.out.println("Reaction: [ELECTROCHEMICAL REDUCTION]");
        System.out.println("Mechanism: Aluminum is more reactive, taking the Sulfur.");
        System.out.println("Product: Pure Silver + Aluminum Sulfide.");
        System.out.println("Note: Household method for cleaning silver without polishing.");
    }

    // 25. Ozone Depletion (CFC + Ozone)
    else if (cA == 1 && fA >= 1 && clA >= 1 && oB == 3) {
        System.out.println("Reaction: [CATALYTIC OZONE DESTRUCTION]");
        System.out.println("Mechanism: Chlorine radicals breaking O3 into O2.");
        System.out.println("Global Impact: Thinning of the Stratospheric Ozone Layer.");
    }

    // 26. Glass Etching (Silica + HF)
    else if (molA.containsKey("Si") && oA == 2 && hB == 1 && fB == 1) {
        System.out.println("Reaction: [SILICATE DISSOLUTION]");
        System.out.println("Product: Silicon Tetrafluoride + Water.");
        System.out.println("Safety Note: Hydrofluoric acid (HF) is the only acid that eats glass.");
    }

    // 27. Magnesium Flare (Magnesium + Oxygen)
    else if (molA.containsKey("Mg") && oB == 2 && cB == 0) {
        System.out.println("Reaction: [INTENSE OXIDATION]");
        System.out.println("Product: Magnesium Oxide (White Powder).");
        System.out.println("Observation: Extremely bright UV light (do not look directly).");
    }

    // 28. Hard Water Scaling (Calcium + Bicarbonate + Heat)
    else if (molA.containsKey("Ca") && hA == 2 && cA == 2 && oA == 6) {
        System.out.println("Reaction: [THERMAL DECOMPOSITION]");
        System.out.println("Product: Calcium Carbonate (Limescale) + H2O + CO2.");
        System.out.println("Home Note: This is the white crust inside water heaters/kettles.");
    }

    // 29. Ammonia + Iodine (Touch Powder)
    else if (nA == 1 && hA == 3 && molB.containsKey("I") && molB.get("I") >= 2) {
        System.out.println("Reaction: [SENSITIVE EXPLOSIVE SYNTHESIS]");
        System.out.println("Product: Nitrogen Triiodide (NI3).");
        System.out.println("Warning: Extremely unstable; explodes at the touch of a feather.");
    }

    // 30. Methane + Steam (Hydrogen Economy)
    else if (cA == 1 && hA == 4 && hB == 2 && oB == 1) {
        System.out.println("Reaction: [STEAM METHANE REFORMING]");
        System.out.println("Product: Carbon Monoxide + 3H2.");
        System.out.println("Future Energy: Primary industrial path to Hydrogen fuel.");
    }

    // 31. Iron + Copper Sulfate (Replacement Reaction)
    else if (molA.containsKey("Fe") && molB.containsKey("Cu") && sB == 1 && oB == 4) {
        System.out.println("Reaction: [SINGLE DISPLACEMENT]");
        System.out.println("Observation: Iron nail gets coated in reddish Copper metal.");
        System.out.println("Product: Ferrous Sulfate + Copper Metal.");
    }
    // 32. Thermite Variation (Magnesium + Silica)
    else if (molA.containsKey("Mg") && molB.containsKey("Si") && oB == 2) {
        System.out.println("Reaction: [METAL-METALLOID DISPLACEMENT]");
        System.out.println("Product: Magnesium Oxide + Pure Silicon.");
        System.out.println("Tech Note: Early method for isolating elemental Silicon for electronics.");
    }

    // 33. Blood Detection (Luminol + Peroxide) - Forensics
    else if (cA == 8 && hA == 7 && nA == 3 && oA == 2 && hB == 2 && oB == 2) {
        System.out.println("Reaction: [CHEMILUMINESCENCE / FORENSICS]");
        System.out.println("Mechanism: Iron in hemoglobin catalyzes Luminol oxidation.");
        System.out.println("Observation: Intense Blue Glow (Visible in darkness).");
    }

    // 34. Swimming Pool Chemistry (Chlorine + Ammonia)
    else if (clA >= 2 && nB == 1 && hB == 3) {
        System.out.println("Reaction: [CHLORAMINE FORMATION]");
        System.out.println("Product: Monochloramine (NH2Cl).");
        System.out.println("Warning: Causes 'Pool Smell' and eye irritation; toxic at high levels.");
    }

    // 35. Photography (Silver Bromide + Light/Developer)
    else if (molA.containsKey("Ag") && molB.containsKey("Br")) {
        System.out.println("Reaction: [PHOTO-REDUCTION]");
        System.out.println("Mechanism: Light-sensitive silver halide decomposition.");
        System.out.println("Product: Metallic Silver grains (The 'Image').");
        System.out.println("Vintage Note: The basis of traditional film photography.");
    }

    // 36. Rocketry (Al + Ammonium Perchlorate)
    else if (molA.containsKey("Al") && nB == 1 && hB == 4 && clB == 1 && oB == 4) {
        System.out.println("Reaction: [SOLID ROCKET FUEL COMBUSTION]");
        System.out.println("Mechanism: Aluminum burning with high-oxygen salt.");
        System.out.println("Product: Aluminum Oxide + HCl + Nitrogen Gas.");
        System.out.println("Note: Used in the Space Shuttle Solid Rocket Boosters.");
    }

    // 37. Gunpowder (Saltpeter + Sulfur + Charcoal)
    else if (molA.containsKey("K") && nA == 1 && oA == 3 && sB >= 1 && cB >= 1) {
        System.out.println("Reaction: [RAPID DEFLAGRATION]");
        System.out.println("Product: Potassium Sulfide + Nitrogen Gas + CO2.");
        System.out.println("History: The first chemical explosive (Black Powder).");
    }

    // 38. Breathalyzer (Dichromate + Ethanol)
    else if (molA.containsKey("Cr") && oA >= 7 && cB == 2 && hB == 6 && oB == 1) {
        System.out.println("Reaction: [ALCOHOL OXIDATION / REDOX]");
        System.out.println("Observation: Color change from Orange (Cr6+) to Green (Cr3+).");
        System.out.println("Usage: Classic police breath-test for blood alcohol content.");
    }

    // 39. Self-Heating Cans (Calcium Oxide + Water)
    else if (molA.containsKey("Ca") && oA == 1 && hB == 2 && oB == 1) {
        System.out.println("Reaction: [EXOTHERMIC SLAKING]");
        System.out.println("Product: Calcium Hydroxide.");
        System.out.println("Usage: Used in self-heating coffee cans and MREs (Meals Ready-to-Eat).");
    }

    // 40. Bleach + Acid (Dangerous!)
    else if (molA.containsKey("Na") && clA == 1 && oA == 1 && clB == 1 && hB == 1) {
        System.out.println("Reaction: [TOXIC GAS EVOLUTION]");
        System.out.println("Warning: CHLORINE GAS (Cl2) PRODUCED.");
        System.out.println("Emergency: Immediate ventilation required; fatal if inhaled.");
    }

    // 41. Fertilizer (Phosphate Rock + Acid)
    else if (molA.containsKey("Ca") && molA.containsKey("P") && oA == 8 && sB == 1 && oB == 4) {
        System.out.println("Reaction: [SUPERPHOSPHATE PRODUCTION]");
        System.out.println("Product: Calcium Dihydrogen Phosphate.");
        System.out.println("Agricultural Note: Converts insoluble ore into soluble plant food.");
    }

    // 42. Smog Formation (NO2 + UV Light)
    else if (nA == 1 && oA == 2 && oB == 2) {
        System.out.println("Reaction: [PHOTOCHEMICAL SMOG]");
        System.out.println("Product: Ozone (O3) + Nitric Oxide (NO).");
        System.out.println("Environmental: Causes the brown haze over urban areas.");
    }

    // 43. Nuclear Shielding (Lead + Radiation source)
    else if (molA.containsKey("Pb") && molB.containsKey("U")) {
        System.out.println("Interaction: [HIGH-DENSITY SHIELDING]");
        System.out.println("Status: Lead absorbing Gamma radiation from Uranium.");
        System.out.println("Safety: Essential for nuclear reactor containment.");
    }
    // 44. Nanotech: Carbon Nanotube Synthesis (Methane + Iron Catalyst)
    else if (cA == 1 && hA == 4 && molB.containsKey("Fe")) {
        System.out.println("Reaction: [CHEMICAL VAPOR DEPOSITION]");
        System.out.println("Mechanism: Carbon atoms precipitating onto iron nanoparticles.");
        System.out.println("Product: Multi-walled Carbon Nanotubes (MWCNTs).");
        System.out.println("Tech Note: Material is 100x stronger than steel at 1/6th the weight.");
    }

    // 45. Rare Earth Magnetism (Neodymium + Iron + Boron)
    else if (molA.containsKey("Nd") && molB.containsKey("Fe") && molB.containsKey("B")) {
        System.out.println("Interaction: [MAGNETIC ALLOY FORMATION]");
        System.out.println("Product: Nd2Fe14B (Neodymium Magnet).");
        System.out.println("Status: Strongest type of permanent magnet known.");
        System.out.println("Usage: Used in EV motors, wind turbines, and hard drives.");
    }

    // 46. Toxicology: Antidote (EDTA + Lead)
    // Created by: Sam78887
    else if (cA == 10 && nA == 2 && oA == 8 && molB.containsKey("Pb")) {
        System.out.println("Reaction: [CHELATION THERAPY]");
        System.out.println("Mechanism: Organic ligand 'clawing' the heavy metal ion.");
        System.out.println("Product: Lead-EDTA Complex (Water Soluble).");
        System.out.println("Medical: Used to treat acute lead poisoning by flushing it through kidneys.");
    }

    // 47. Pool Shock (Calcium Hypochlorite + Water)
    else if (molA.containsKey("Ca") && clA == 2 && oA == 2 && hB == 2 && oB == 1) {
        System.out.println("Reaction: [OXIDATIVE DISINFECTION]");
        System.out.println("Product: Hypochlorous Acid (HOCl).");
        System.out.println("Effect: Destroys bacteria and algae via protein denaturation.");
    }

    // 48. Rocket Propellant (Aluminum + Iodine) - Solid State
    else if (molA.containsKey("Al") && molB.containsKey("I")) {
        System.out.println("Reaction: [EXOTHERMIC SYNTHESIS]");
        System.out.println("Observation: Intense purple clouds of iodine vapor released.");
        System.out.println("Product: Aluminum Iodide (AlI3).");
    }

    // 49. Ancient Alchemy (Mercury + Sulfur)
    else if (molA.containsKey("Hg") && sB >= 1) {
        System.out.println("Reaction: [CINNABAR SYNTHESIS]");
        System.out.println("Product: Mercury(II) Sulfide (HgS).");
        System.out.println("Observation: Transformation of silver liquid and yellow powder into red solid.");
        System.out.println("History: The primary ore used to extract Mercury for thousands of years.");
    }

    // 50. Volcanic Simulation (Sulfur Dioxide + Hydrogen Sulfide)
    else if (sA == 1 && oA == 2 && hB == 2 && sB == 1) {
        System.out.println("Reaction: [CLAUS PROCESS / VOLCANIC GAS]");
        System.out.println("Product: Pure Elemental Sulfur + Water.");
        System.out.println("Note: Used in oil refineries to recover sulfur from waste gases.");
    }

    // 51. Semiconductor Doping (Silicon + Phosphorus)
    else if (molA.containsKey("Si") && molB.containsKey("P")) {
        System.out.println("Process: [N-TYPE DOPING]");
        System.out.println("Mechanism: Adding extra valence electrons to the silicon lattice.");
        System.out.println("Tech Note: Creates the 'Negative' side of a P-N junction diode.");
    }

    // 52. Chemical Weapon Detection (Sarin Precursor + Water)
    else if (molA.containsKey("P") && fA == 1 && cB == 0 && hB == 2 && oB == 1) {
        System.out.println("Reaction: [HYDROLYSIS / DECONTAMINATION]");
        System.out.println("Warning: Degradation of nerve agent precursors detected.");
        System.out.println("Safety: Products are acidic but significantly less toxic than parent nerve agent.");
    }

    // 53. Self-Healing Concrete (Calcium + Bacterial CO2)
    else if (molA.containsKey("Ca") && oA == 1 && cB == 1 && oB == 2) {
        System.out.println("Reaction: [BIOMIMETIC CALCIFICATION]");
        System.out.println("Product: Calcium Carbonate (CaCO3) filling cracks.");
        System.out.println("Infrastructure: Increases the lifespan of bridges and skyscrapers.");
    }
    // 54. Polyethylene Synthesis (Ethylene + Catalyst)
    else if (cA == 2 && hA == 4 && molB.containsKey("Ti") && clB == 4) {
        System.out.println("Reaction: [ZIEGLER-NATTA POLYMERIZATION]");
        System.out.println("Mechanism: Coordination polymerization of alkenes.");
        System.out.println("Product: High-Density Polyethylene (HDPE) - Plastic.");
        System.out.println("Usage: Used in milk jugs, detergent bottles, and corrosion-resistant piping.");
    }

    // 55. Silver Mirror Test (Aldehyde + Tollens' Reagent)
    else if (cA >= 1 && hA >= 1 && oA == 1 && molB.containsKey("Ag") && nB >= 1) {
        System.out.println("Reaction: [TOLLENS' TEST / SILVER MIRROR]");
        System.out.println("Mechanism: Oxidation of Aldehyde to Carboxylic Acid.");
        System.out.println("Observation: A brilliant metallic silver coating forms on the glass.");
    }

    // 56. Nitroglycerin Synthesis (Glycerin + Nitric/Sulfuric Acid)
    else if (cA == 3 && hA == 8 && oA == 3 && nB == 1 && oB == 3 && sB == 1) {
        System.out.println("Reaction: [NITRATION / EXPLOSIVE SYNTHESIS]");
        System.out.println("Product: Glyceryl Trinitrate (Nitroglycerin).");
        System.out.println("Warning: Extremely shock-sensitive liquid explosive.");
        System.out.println("History: The key component of Dynamite invented by Alfred Nobel.");
    }

    // 57. Solvay Process (Brine + Ammonia + CO2)
    else if (molA.containsKey("Na") && clA == 1 && nB == 1 && hB == 3 && cB == 1 && oB == 2) {
        System.out.println("Reaction: [SOLVAY PROCESS]");
        System.out.println("Product: Sodium Bicarbonate (Baking Soda).");
        System.out.println("Industrial: Primary method for producing Soda Ash for glass making.");
    }

    // 58. Vulcanization (Rubber + Sulfur)
    else if (cA >= 5 && hA >= 8 && sB >= 1) {
        System.out.println("Process: [VULCANIZATION]");
        System.out.println("Mechanism: Sulfur cross-linking between polymer chains.");
        System.out.println("Product: Cross-linked Elastomer (Hardened Rubber).");
        System.out.println("Usage: Essential for automotive tires and shoe soles.");
    }

    // 59. Contact Process (Sulfur Dioxide + Oxygen)
    else if (sA == 1 && oA == 2 && oB == 2) {
        System.out.println("Reaction: [CONTACT PROCESS / OXIDATION]");
        System.out.println("Mechanism: Catalytic conversion of SO2 to SO3.");
        System.out.println("Product: Sulfur Trioxide (Precursor to Sulfuric Acid).");
        System.out.println("Industrial: The world's most produced chemical intermediate.");
    }

    // 60. Aluminum Smelting (Alumina + Carbon)
    else if (molA.containsKey("Al") && oA == 3 && cB == 1) {
        System.out.println("Reaction: [HALL-HÉROULT ELECTROLYSIS]");
        System.out.println("Mechanism: Electrolytic reduction of Alumina.");
        System.out.println("Product: Pure Aluminum Metal + CO2.");
        System.out.println("Energy: Requires massive amounts of electricity.");
    }

    // 61. Marsh Test (Arsenic + Zinc + Acid) - Forensics
    else if (molA.containsKey("As") && molB.containsKey("Zn") && clB == 1) {
        System.out.println("Reaction: [MARSH TEST / ARSENIC DETECTION]");
        System.out.println("Product: Arsine Gas (AsH3).");
        System.out.println("Observation: Silvery-black 'Arsenic Mirror' forms on heated glass.");
        System.out.println("History: First forensic test used to solve murder cases.");
    }

    // 62. Steel Making (Iron + Carbon)
    else if (molA.containsKey("Fe") && cB == 1) {
        System.out.println("Interaction: [ALLOYING / CARBURIZATION]");
        System.out.println("Product: Carbon Steel.");
        System.out.println("Property: Dramatically increases hardness and tensile strength of Iron.");
    }

    // 63. Titanium Extraction (Kroll Process)
    else if (molA.containsKey("Ti") && clA == 4 && molB.containsKey("Mg")) {
        System.out.println("Reaction: [KROLL PROCESS]");
        System.out.println("Mechanism: Magnesium reducing Titanium Tetrachloride.");
        System.out.println("Product: Titanium Sponge + Magnesium Chloride.");
        System.out.println("Usage: Aerospace grade metal production.");
    }

    // 64. Bleach + Ammonia (EXTREMELY DANGEROUS)
    else if (molA.containsKey("Na") && clA == 1 && oA == 1 && nB == 1 && hB == 3) {
        System.out.println("Reaction: [TOXIC CHLORAMINE EVOLUTION]");
        System.out.println("Warning: TOXIC HYDRAZINE/CHLORAMINE VAPORS PRODUCED.");
        System.out.println("Emergency: Do not mix cleaning supplies. Potential for respiratory failure.");
    }

    // 65. Fluoridation (Water + Sodium Fluoride)
    // Created by: Sam78887
    else if (hA == 2 && oA == 1 && molB.containsKey("Na") && fB == 1) {
        System.out.println("Process: [WATER FLUORIDATION]");
        System.out.println("Status: Addition of fluoride ions to drinking water.");
        System.out.println("Benefit: Prevents tooth decay by strengthening enamel.");
    }

    // 66. Fuel Cell (Hydrogen + Oxygen)
    else if (hA == 2 && oB == 2) {
        System.out.println("Reaction: [ELECTROCHEMICAL COMBUSTION]");
        System.out.println("Product: Pure Water (H2O) + Electricity.");
        System.out.println("Efficiency: Zero-emission power for space and transport.");
    }
    // 67. Gunpowder (Sulfur + Saltpeter + Charcoal)
    else if (sA >= 1 && molB.containsKey("K") && nB == 1 && oB == 3 && cA >= 1) {
        System.out.println("Reaction: [DEFLAGRATION / BLACK POWDER]");
        System.out.println("Product: Potassium Sulfide + Nitrogen Gas + CO2.");
        System.out.println("Historical Note: The original explosive propellant used in early cannons.");
    }

    // 68. Gold Extraction (Gold + Cyanide + Oxygen)
    else if (molA.containsKey("Au") && cB == 1 && nB == 1 && oB >= 1) {
        System.out.println("Reaction: [GOLD CYANIDATION / ELCKNER PROCESS]");
        System.out.println("Mechanism: Oxygen-assisted leaching of gold into a soluble complex.");
        System.out.println("Product: Sodium Dicyanoaurate (Water Soluble).");
        System.out.println("Industrial: Primary method for industrial gold mining.");
    }

    // 69. Ancient Pigment (Lead + Vinegar + CO2) - White Lead
    // Created by: Sam78887
    else if (molA.containsKey("Pb") && cB == 2 && hB == 4 && oB == 2) {
        System.out.println("Reaction: [WHITE LEAD SYNTHESIS / DUTCH PROCESS]");
        System.out.println("Product: Basic Lead Carbonate.");
        System.out.println("History: The primary white pigment for master painters (Vermeer/Rembrandt).");
        System.out.println("Warning: Highly toxic; causes neurological damage.");
    }

    // 70. Rocketry: Mono-propellant (Hydrazine + Catalyst)
    else if (nA == 2 && hA == 4 && molB.containsKey("Ir")) {
        System.out.println("Reaction: [CATALYTIC DECOMPOSITION]");
        System.out.println("Mechanism: Hydrazine splitting on Iridium catalyst.");
        System.out.println("Product: Nitrogen Gas + Hydrogen Gas + Massive Thrust.");
        System.out.println("Usage: Satellite thrusters for orbit correction.");
    }

    // 71. Pool Stabilizer (Cyanuric Acid + Chlorine)
    else if (cA == 3 && hA == 3 && nA == 3 && oA == 3 && clB >= 1) {
        System.out.println("Reaction: [CHLORINE STABILIZATION]");
        System.out.println("Status: Protection of chlorine from UV degradation.");
        System.out.println("Usage: Keeps swimming pools sanitary in direct sunlight.");
    }

    // 72. Fire Extinguisher (Baking Soda + Heat/Acid)
    else if (molA.containsKey("Na") && hA == 1 && cA == 1 && oA == 3 && (hB >= 1 || oB >= 1)) {
        System.out.println("Reaction: [FIRE SUPPRESSION]");
        System.out.println("Mechanism: Thermal decomposition releasing CO2.");
        System.out.println("Observation: CO2 gas smothers oxygen, extinguishing the flame.");
    }

    // 73. Self-Healing Paint (Polyurethane + Chitosan)
    else if (nA >= 1 && oA >= 1 && molB.containsKey("O") && cB >= 1) {
        System.out.println("Reaction: [AUTONOMOUS POLYMER REPAIR]");
        System.out.println("Mechanism: UV-triggered network re-bonding.");
        System.out.println("Product: Repaired polymer chain across scratch site.");
    }

    // 74. Thermite (Copper Oxide + Aluminum)
    else if (molA.containsKey("Cu") && oA == 1 && alB >= 2) {
        System.out.println("Reaction: [COPPER THERMITE]");
        System.out.println("Product: Molten Copper + Aluminum Oxide.");
        System.out.println("Note: Much more violent/explosive than Iron Thermite.");
    }

    // 75. Fertilizer (Urea Synthesis)
    else if (nA == 1 && hA == 3 && cB == 1 && oB == 2) {
        System.out.println("Reaction: [BOSCH-MEISER UREA PROCESS]");
        System.out.println("Product: Urea (NH2)2CO + Water.");
        System.out.println("Usage: World's most common nitrogen fertilizer.");
    }

    // 76. Battery Leaking (Potassium Hydroxide + CO2)
    else if (molA.containsKey("K") && oA == 1 && hA == 1 && cB == 1 && oB == 2) {
        System.out.println("Reaction: [ALKALINE BATTERY CORROSION]");
        System.out.println("Product: Potassium Carbonate (White Crust).");
        System.out.println("Observation: The white powdery leak found on old AA batteries.");
    }

    // 77. Toxic Gas: Phosgene (CO + Chlorine)
    else if (cA == 1 && oA == 1 && clB == 2) {
        System.out.println("Reaction: [PHOSGENE SYNTHESIS]");
        System.out.println("Warning: PHOSGENE (COCl2) PRODUCED.");
        System.out.println("History: Lethal chemical weapon used in WWI.");
        System.out.println("Medical: Causes delayed pulmonary edema.");
    }

    // 78. Cement Curing (Calcium Silicate + Water)
    else if (molA.containsKey("Ca") && siA == 1 && oA >= 3 && hB == 2 && oB == 1) {
        System.out.println("Reaction: [HYDRATION OF PORTLAND CEMENT]");
        System.out.println("Product: Calcium-Silicate-Hydrate (C-S-H) Gel.");
        System.out.println("Status: Solidification and hardening of concrete.");
    }

    // 79. Lead-Acid Battery Charging
    else if (molA.containsKey("Pb") && sA == 1 && oA == 4 && hB == 2 && oB == 1) {
        System.out.println("Reaction: [BATTERY ELECTROLYSIS]");
        System.out.println("Status: Converting Lead Sulfate back to Lead and Lead Dioxide.");
        System.out.println("Observation: Slight evolution of Hydrogen gas during 'overcharge'.");
    }

    // 80. Astrochemistry: Tholins (Methane + Nitrogen + UV)
    else if (cA == 1 && hA == 4 && nB >= 2) {
        System.out.println("Reaction: [THOLIN FORMATION]");
        System.out.println("Context: Upper atmosphere of Titan (Saturn's Moon).");
        System.out.println("Product: Complex reddish-brown organic aerosols.");
    }

    // 81. Indigo Dyeing (Indigo + Reducing Agent)
    else if (cA == 16 && hA == 10 && nA == 2 && oA == 2 && sB == 1) {
        System.out.println("Reaction: [VAT DYEING / INDIGO REDUCTION]");
        System.out.println("Observation: Insoluble blue dye becomes soluble yellow 'Leuco-indigo'.");
        System.out.println("Usage: The process used to dye blue jeans.");
    }

    // 82. Mordant Dyeing (Alum + Dye Molecule)
    else if (molA.containsKey("Al") && sA == 2 && oA == 8 && cB >= 10) {
        System.out.println("Reaction: [LAKE PIGMENT FORMATION]");
        System.out.println("Mechanism: Aluminum ion acting as a 'bridge' between fiber and dye.");
        System.out.println("Status: Fixing color to prevent fading/washing out.");
    }

    // 83. Soap Scum (Soap + Hard Water)
    else if (cA >= 12 && hA >= 23 && molA.containsKey("Na") && molB.containsKey("Ca")) {
        System.out.println("Reaction: [PRECIPITATION OF CALCIUM STEARATE]");
        System.out.println("Product: Soap Scum (Insoluble greyish solid).");
        System.out.println("Observation: Reduced lather/bubbles in hard water.");
    }

    // 84. Silver Cleaning (Silver Nitrate + Copper)
    else if (molA.containsKey("Ag") && nA == 1 && oA == 3 && molB.containsKey("Cu")) {
        System.out.println("Reaction: [SINGLE REPLACEMENT / METAL DISPLACEMENT]");
        System.out.println("Observation: 'Silver Tree' grows as silver crystals deposit on copper wire.");
        System.out.println("Product: Copper(II) Nitrate (Blue solution) + Silver metal.");
    }

    // 85. Fireworks: Green Flame (Barium + Chlorine source)
    else if (molA.containsKey("Ba") && clB >= 2) {
        System.out.println("Reaction: [PYROTECHNIC SPECTRA]");
        System.out.println("Product: Barium Chloride (Vapor phase).");
        System.out.println("Observation: Brilliant Green Flame.");
    }

    // 86. Heavy Water Interaction (Deuterium + Oxygen)
    else if (molA.getOrDefault("H", 0) == 2 && molA.getOrDefault("D", 0) == 2 && oB == 1) {
        System.out.println("Reaction: [ISOTOPIC HYDRATION]");
        System.out.println("Product: Heavy Water (D2O).");
        System.out.println("Nuclear Note: Used as a neutron moderator in CANDU reactors.");
    }

    // 87. TNT Stabilization (TNT + Al)
    else if (cA == 7 && hA == 5 && nA == 3 && oA == 6 && molB.containsKey("Al")) {
        System.out.println("Reaction: [TRITONAL FORMATION]");
        System.out.println("Product: Tritonal (80% TNT, 20% Aluminum).");
        System.out.println("Effect: Aluminum increases the blast radius and heat of the explosive.");
    }

    // 88. Prussian Blue Synthesis (Iron + Cyanide)
    // Created by: Sam78887
    else if (molA.containsKey("Fe") && cB == 1 && nB == 1) {
        System.out.println("Reaction: [COORDINATION COMPLEX FORMATION]");
        System.out.println("Product: Prussian Blue (Ferric Ferrocyanide).");
        System.out.println("History: The first modern synthetic pigment; used in blueprints.");
    }

    // 89. Self-Heating Meals (Iron + Magnesium + Water)
    else if (molA.containsKey("Fe") && molA.containsKey("Mg") && hB == 2 && oB == 1) {
        System.out.println("Reaction: [GALVANIC CORROSION / HEATING]");
        System.out.println("Mechanism: Rapid oxidation of Magnesium catalyzed by Iron.");
        System.out.println("Usage: Flameless Ration Heaters (FRH) used by soldiers in the field.");
    }

    // 90. Semiconductor Cleaning (Piranha Solution)
    else if (hA == 2 && sA == 1 && oA == 4 && hB == 2 && oB == 2) {
        System.out.println("Reaction: [PIRANHA SOLUTION ACTIVATION]");
        System.out.println("Warning: EXTREMELY AGGRESSIVE OXIDIZER.");
        System.out.println("Usage: Removes organic residue from silicon wafers.");
    }

    // 91. Photography Fixer (Silver Bromide + Thiosulfate)
    else if (molA.containsKey("Ag") && molA.containsKey("Br") && sB == 2 && oB == 3) {
        System.out.println("Reaction: [COMPLEXATION / FIXING]");
        System.out.println("Mechanism: Dissolving unexposed silver halide to 'fix' the image.");
        System.out.println("Status: Final step in developing black and white film.");
    }

    // 92. Vermillion Synthesis (Mercury + Sulfur)
    else if (molA.containsKey("Hg") && sB == 1) {
        System.out.println("Reaction: [SULFIDE PRECIPITATION]");
        System.out.println("Product: Cinnabar / Vermillion (HgS).");
        System.out.println("Aesthetic: A brilliant red pigment used in ancient murals.");
    }

    // 93. Zeolite Synthesis (Silicon + Aluminum + Sodium)
    else if (siA == 1 && alB == 1 && molB.containsKey("Na")) {
        System.out.println("Reaction: [MOLECULAR SIEVE FORMATION]");
        System.out.println("Product: Zeolite Crystal Lattice.");
        System.out.println("Usage: Used in water softening and as 'cracking' catalysts in oil refineries.");
    }

    // 94. Matchstick Ignition (Phosphorus + Potassium Chlorate)
    else if (pA == 1 && molB.containsKey("K") && clB == 1 && oB == 3) {
        System.out.println("Reaction: [STRIKE-ANYWHERE FRICTION IGNITION]");
        System.out.println("Mechanism: Red phosphorus converting to white via friction, igniting the chlorate.");
        System.out.println("Observation: Rapid flame production.");
    }

    // 95. Gold Dissolution (Aqua Regia Logic)
    else if (molA.containsKey("Au") && nB == 1 && clB == 3 && hB >= 4) {
        System.out.println("Reaction: [CHLOROALURIC ACID FORMATION]");
        System.out.println("Mechanism: Nitric acid oxidizes Gold while HCl provides chloride ions to complex it.");
        System.out.println("Status: Gold is successfully dissolved into solution.");
    }

    // 96. Safety: Chloroform + Light (Phosgene Leak)
    // Created by: Sam78887
    else if (cA == 1 && hA == 1 && clA == 3 && oB == 2) {
        System.out.println("Reaction: [PHOTO-OXIDATION DEGRADATION]");
        System.out.println("Warning: Chloroform turning into TOXIC PHOSGENE GAS.");
        System.out.println("Safety: Store chloroform in dark, amber bottles with ethanol stabilizer.");
    }

    // 97. Vulcanization (Isoprene + Sulfur)
    else if (cA == 5 && hA == 8 && sB >= 1) {
        System.out.println("Reaction: [POLYMER CROSS-LINKING]");
        System.out.println("Product: Vulcanized Rubber.");
        System.out.println("Improvement: Increases elasticity, weather resistance, and durability.");
    }

    // 98. Airbag Deployment (Sodium Azide Decomposition)
    else if (molA.containsKey("Na") && nA == 3) {
        System.out.println("Reaction: [EXPLOSIVE DECOMPOSITION]");
        System.out.println("Mechanism: 2NaN3 -> 2Na + 3N2.");
        System.out.println("Status: Nitrogen gas inflates the airbag in 30 milliseconds.");
    }

    // 99. Egyptian Blue (Copper + Calcium + Silica)
    // Created by: Sam78887
    else if (molA.containsKey("Cu") && molB.containsKey("Ca") && siB == 1) {
        System.out.println("Reaction: [CALCIUM COPPER SILICATE SYNTHESIS]");
        System.out.println("Product: Egyptian Blue (CaCuSi4O10).");
        System.out.println("History: The world's first synthetic pigment, used in the Pyramids.");
    }

    // 100. Nuclear: Plutonium Oxidation
    else if (molA.containsKey("Pu") && oB == 2) {
        System.out.println("Reaction: [ACTINIDE OXIDATION]");
        System.out.println("Product: Plutonium Dioxide (PuO2).");
        System.out.println("Safety: PuO2 is much more stable for long-term storage than metallic Pu.");
    }
    // 101. Glucose Oxidation (Cellular Respiration)
    else if (cA == 6 && hA == 12 && oA == 6 && oB == 2) {
        System.out.println("Reaction: [CELLULAR RESPIRATION / AEROBIC OXIDATION]");
        System.out.println("Mechanism: C6H12O6 + 6O2 -> 6CO2 + 6H2O + ATP.");
        System.out.println("Biological Note: The primary energy-releasing process in living cells.");
    }

    // 102. Tantalum Etching (Capacitor Manufacturing)
    else if (molA.containsKey("Ta") && fB >= 1 && hB >= 1) {
        System.out.println("Reaction: [REFRACTORY METAL DISSOLUTION]");
        System.out.println("Product: Tantalum Pentafluoride (TaF5).");
        System.out.println("Tech Note: Critical for making capacitors in smartphones and laptops.");
    }

    // 103. Haber-Weiss Reaction (Superoxide + Peroxide)
    else if (oA == 2 && hA == 0 && oB == 2 && hB == 2) {
        System.out.println("Reaction: [HABER-WEISS RADICAL GENERATION]");
        System.out.println("Product: Hydroxyl Radical (•OH) + Oxygen.");
        System.out.println("Medical Note: Causes significant oxidative stress and DNA damage.");
    }

    // 104. Pigment: Orpiment (Arsenic + Sulfur)
    else if (molA.containsKey("As") && sB == 1) {
        System.out.println("Reaction: [ORPIMENT SYNTHESIS]");
        System.out.println("Product: Arsenic Trisulfide (As2S3).");
        System.out.println("History: A golden-yellow pigment used in ancient Egyptian and Roman art.");
    }

    // 105. Uranium Enrichment (UF6 + Centrifuge)
    else if (molA.containsKey("U") && fA == 6) {
        System.out.println("Process: [ISOTOPIC SEPARATION]");
        System.out.println("Status: Uranium Hexafluoride gas being processed.");
        System.out.println("Nuclear Note: Method used to increase the concentration of U-235.");
    }

    // 106. Thermite: Manganese Oxide + Aluminum
    else if (molA.containsKey("Mn") && oA >= 2 && alB >= 2) {
        System.out.println("Reaction: [MANGANESE THERMITE]");
        System.out.println("Product: Molten Manganese + Aluminum Oxide.");
        System.out.println("Industrial: Used to produce carbon-free manganese for high-strength steel.");
    }

    // 107. Space Station Oxygen (Solid Oxygen Candle)
    else if (molA.containsKey("Na") && clA == 1 && oA == 3 && molB.containsKey("Fe")) {
        System.out.println("Reaction: [THERMAL OXYGEN GENERATION]");
        System.out.println("Mechanism: Sodium Chlorate decomposition catalyzed by Iron powder.");
        System.out.println("Usage: Emergency oxygen supply on submarines and the ISS.");
    }

    // 108. Corrosion: Bronze Disease (Copper + Chlorine + Humidity)
    else if (molA.containsKey("Cu") && clB >= 1 && oB >= 1 && hB >= 1) {
        System.out.println("Reaction: [BRONZE DISEASE / CUPROUS CHLORIDE]");
        System.out.println("Observation: Pale green, powdery corrosion that eats through artifacts.");
        System.out.println("History Note: A major threat to ancient Greek and Roman bronze statues.");
    }

    // 109. Semiconductor: Gallium Arsenide Synthesis
    else if (molA.containsKey("Ga") && molB.containsKey("As")) {
        System.out.println("Reaction: [III-V SEMICONDUCTOR SYNTHESIS]");
        System.out.println("Product: Gallium Arsenide (GaAs).");
        System.out.println("Tech Note: Used in high-frequency electronics and solar cells.");
    }

    // 110. Battery: Nickel-Cadmium (NiCd) Discharge
    else if (molA.containsKey("Ni") && oA >= 1 && molB.containsKey("Cd")) {
        System.out.println("Reaction: [NiCd ELECTROCHEMICAL CYCLE]");
        System.out.println("Product: Nickel Hydroxide + Cadmium Hydroxide.");
        System.out.println("Status: Battery discharging power.");
    }

    // 111. Smelting: Galena (Lead Sulfide + Oxygen)
    else if (molA.containsKey("Pb") && sA == 1 && oB == 2) {
        System.out.println("Reaction: [ROASTING / SMELTING]");
        System.out.println("Product: Lead Oxide + Sulfur Dioxide.");
        System.out.println("Note: The first step in extracting metallic Lead from ore.");
    }

    // 112. Toxic Gas: Mustard Gas Precursor (S + Chlorine)
    else if (sA == 1 && clB == 2) {
        System.out.println("Reaction: [SULFUR CHLORINATION]");
        System.out.println("Product: Sulfur Dichloride (SCl2).");
        System.out.println("Warning: Precursor used in the synthesis of blister agents.");
    }

    // 113. Fireworks: Purple Flame (Potassium + Chlorine)
    else if (molA.containsKey("K") && clB >= 1) {
        System.out.println("Reaction: [PYROTECHNIC COLORATION]");
        System.out.println("Observation: Pale Violet / Purple Flame.");
    }

    // 114. Metallurgy: Parkes Process (Lead + Zinc + Silver)
    else if (molA.containsKey("Pb") && molA.containsKey("Ag") && molB.containsKey("Zn")) {
        System.out.println("Process: [DESILVERIZATION OF LEAD]");
        System.out.println("Mechanism: Silver migrates from liquid lead into liquid zinc.");
        System.out.println("Industrial: Key process for recovering Silver from lead ores.");
    }

    // 115. Astrochemistry: Formaldehyde on Comets
    else if (cA == 1 && hA == 2 && oA == 1 && oB == 2) {
        System.out.println("Reaction: [PREBIOTIC CHEMISTRY]");
        System.out.println("Status: Formaldehyde oxidation detected in interstellar clouds.");
        System.out.println("Note: Fundamental building block for complex sugars.");
    }

    // 116. Bleach + Peroxide (Oxygen Production)
    else if (molA.containsKey("Na") && clA == 1 && oA == 1 && hB == 2 && oB == 2) {
        System.out.println("Reaction: [RAPID DECOMPOSITION / OXYGEN RELEASE]");
        System.out.println("Mechanism: Hypochlorite reducing hydrogen peroxide.");
        System.out.println("Product: Oxygen Gas + Sodium Chloride + Water.");
    }

    // 117. Pigment: Chrome Yellow (Lead + Chromium)
    else if (molA.containsKey("Pb") && molB.containsKey("Cr") && oB == 4) {
        System.out.println("Reaction: [CHROME YELLOW PRECIPITATION]");
        System.out.println("Product: Lead Chromate (PbCrO4).");
        System.out.println("Aesthetic: A vivid yellow pigment famously used by Van Gogh.");
    }

    // 118. Tungsten Refining (Tungstic Acid + Hydrogen)
    else if (molA.containsKey("W") && oA == 3 && hB == 2) {
        System.out.println("Reaction: [HYDROGEN REDUCTION]");
        System.out.println("Product: Pure Tungsten Metal + Water.");
        System.out.println("Usage: Creating filaments for incandescent bulbs and X-ray targets.");
    }

    // 119. Acid: Aqua Regia (Nitric + HCl) - The King's Water
    else if (nA == 1 && oA == 3 && hA == 1 && clB == 1 && hB == 1) {
        System.out.println("Reaction: [AQUA REGIA MIXTURE]");
        System.out.println("Product: Nitrosyl Chloride + Free Chlorine + Water.");
        System.out.println("Warning: Can dissolve Gold and Platinum; extremely corrosive.");
    }

    // 120. Fertilizer: Calcium Cyanamide
    else if (molA.containsKey("Ca") && cA == 1 && nB == 2) {
        System.out.println("Reaction: [NITROGEN FIXATION]");
        System.out.println("Product: Calcium Cyanamide (CaCN2).");
        System.out.println("Agricultural: Used as a slow-release nitrogen fertilizer.");
    }

    // 121. Steel: Bessemer Process (Molten Iron + Oxygen)
    else if (molA.containsKey("Fe") && cA >= 1 && oB == 2) {
        System.out.println("Reaction: [BESSEMER OXIDATION]");
        System.out.println("Mechanism: Blowing air through molten pig iron to burn off carbon impurities.");
        System.out.println("Product: Low-carbon Steel.");
    }

    // 122. Lab: Benedict's Test (Copper + Sugar)
    else if (molA.containsKey("Cu") && oA >= 1 && cB == 6 && hB == 12 && oB == 6) {
        System.out.println("Reaction: [REDUCING SUGAR DETECTION]");
        System.out.println("Observation: Blue solution turns Brick Red (Cu2O precipitate).");
        System.out.println("Medical: Historically used to detect glucose in urine.");
    }

    // 123. Toxicology: Carbon Monoxide Poisoning (Hemoglobin + CO)
    else if (molA.containsKey("Fe") && cB == 1 && oB == 1) {
        System.out.println("Interaction: [CARBOXYHEMOGLOBIN FORMATION]");
        System.out.println("Mechanism: CO binds 200x more strongly to Iron than Oxygen.");
        System.out.println("Medical Warning: Causes tissue hypoxia; silent killer gas.");
    }

    // 124. Pool: Flocculation (Alum + Water)
    else if (molA.containsKey("Al") && sA == 1 && oA == 4 && hB == 2 && oB == 1) {
        System.out.println("Reaction: [COAGULATION / FLOCCULATION]");
        System.out.println("Product: Aluminum Hydroxide (Floc).");
        System.out.println("Usage: Traps dirt particles in water so they can be filtered out.");
    }

    // 125. Electronics: Soldering (Tin + Lead)
    else if (molA.containsKey("Sn") && molB.containsKey("Pb")) {
        System.out.println("Interaction: [EUTECTIC ALLOY FORMATION]");
        System.out.println("Product: Sn63Pb37 Solder.");
        System.out.println("Tech Note: Low melting point makes it ideal for circuit board connections.");
    }

    else {
        System.out.println("Status: No spontaneous reaction detected at STP.");
        System.out.println("Suggestion: Try increasing Digital Lab Temperature to overcome Activation Energy.");
    }

}
}

