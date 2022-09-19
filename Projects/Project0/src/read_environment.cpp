// THIS PROGRAM REQUIRES NO MODIFICATION

#include <iostream>
#include <cstdlib>

int main(int argc, char *argv[])
{
    const char *uname = std::getenv("USER");
    const char *gcc = std::getenv("GCC_ROOT");
    const char *nodes = std::getenv("SLURM_JOB_NODELIST");

    if (uname && gcc && nodes && argc == 2)
    {
        std::cout << "Your name is:           " << argv[1] << '\n';
        std::cout << "Your username is:       " << uname << '\n';
        std::cout << "Your gcc root is:       " << gcc << '\n';
        std::cout << "You're connected to:    " << nodes << '\n';
    }
    else
    {
        std::cout << "Something went wrong...\n";
    }
}
