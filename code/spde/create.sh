set -e

### Building the C++ code

cd ~/mres/source/program/build

echo ""
echo "Run the Matlab symbolic code?"
select yn in "Yes" "No"; do
    case $yn in
        Yes ) sh buildvec.sh; break ;;
        No ) break ;;
    esac
done

cd ..
echo ""
echo "Build C++ code?"
select yn in "Yes" "No"; do
    case $yn in
        Yes ) cat Problem_init.c build/built.cpp > Problem.cpp; break ;;
        No ) break ;;
    esac
done

echo ""
echo "Compile C++ code?"
select yn in "Yes" "No"; do 
    case $yn in
        Yes ) make; break ;;
        No ) break ;;
    esac
done

echo ""
echo "***"
