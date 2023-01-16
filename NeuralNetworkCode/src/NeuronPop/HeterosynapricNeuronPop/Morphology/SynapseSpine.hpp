#include <valarray>




class SynapseSpine { //This class is poorly done, yes, but it is my effort to make the legacy part of the code to work with branched code. This is partly self-inflicted, but inheritance is not a good option in this case.
//An alternative option would be inheritance with a bool-check to do a static_cast between the parent and child classes in new code.
//Inheritance was avoided because the required virtualization for pointer usage would be a lot of thankless work
    protected:
    //Legacy variables
    unsigned long preNeuronId{};
    unsigned long postNeuronId{};
    double distToSoma{}; //Use this to build a generalizable distance matrix. In the future this is useful for actual soma distance (computed in a function)
    double lastSpike{};
    double theta{}; // heterosynaptic cooperativity
    double weight{}; //Is this even compatible with negative weights? Yes, as the negative weight comes from J, weight is just a factor to multiply
    unsigned long idInMorpho{}; // id for synapse within its population
    unsigned long idInHCS{}; // id for synapse within its a synapse collection
    //Branched variables
    bool IsBranched{false};
    int branchId{}; //This has to be discrete
    int branchPositionId{}; //This has to be discrete
    int distanceFromNode{}; //This will probably be discrete too, as the distance is id*gap, and gap will be 1 um.
    //int uniqueTreeId{}; Useful only in the indexing of the triangular matrix if implemented

    public:
    //Constructors
    SynapseSpine()=default;
    SynapseSpine(double distToSoma, double lastSpike, double weight);
    SynapseSpine(int distanceFromNode, double lastSpike, double weight, int branchId, int branchPositionId);
    //Methods
    //Bool checks
    bool IsBranchedBool() {return IsBranched;}
    bool SetBranchedTrue() {IsBranched=true;}
    //Getters
    //Legacy
    unsigned long getPreNeuronId(){return preNeuronId;};
    unsigned long getPostNeuronId(){return postNeuronId;};
    double getDistToSoma(){return distToSoma;};
    double getLastSpike(){return lastSpike;};
    double getTheta(){return theta;};
    double getWeight(){return weight;};
    unsigned long getIdInMorpho(){return idInMorpho;};
    unsigned long getIdInHCS(){return idInHCS;};
    //Branched variables
    int getBranchId(){return branchId;};
    int getBranchPositionId(){return branchPositionId;};
    int getDistanceFromNode(){return distanceFromNode;};

    //Setters
    //Legacy
    void setPreNeuronId(unsigned long neuronId){preNeuronId=neuronId;};
    void setPostNeuronId(unsigned long neuronId){postNeuronId=neuronId;};
    void setDistToSoma(double dist){distToSoma=dist;};
    void setLastSpike(double time){lastSpike=time;};
    void setTheta(double thetaIn){theta=thetaIn;};
    void setWeight(double weightIn){weight=weightIn;};
    void setIdInMorpho(unsigned long idIn){idInMorpho=idIn;};
    void setIdInHCS(unsigned long idIn){idInHCS=idIn;};
    //Branched variables
    void setBranchId(int idIn){branchId=idIn;};
    void setBranchPositionId(int positionId){branchPositionId=positionId;};
    void setDistanceFromNode(int distance){distanceFromNode=distance;};


    //Recorder functions
    std::valarray<double> getIndividualSynapticProfileLegacy() const;
    std::valarray<double> getIndividualSynapticProfileBranched() const;
    virtual std::valarray<double> getIndividualSynapticProfile() const;
};